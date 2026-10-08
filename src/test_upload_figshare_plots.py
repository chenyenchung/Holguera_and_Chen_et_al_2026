import hashlib
import sqlite3
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from upload_figshare_plots import Client, parallel_transfers, transfer_file, verify_record


class UploadTests(unittest.TestCase):
    def test_processing_file_does_not_block_another_transfer(self):
        with tempfile.TemporaryDirectory() as directory:
            staging = Path(directory)
            db = sqlite3.connect(staging / "upload.sqlite")
            db.execute("CREATE TABLE files (path TEXT PRIMARY KEY, id INTEGER, md5 TEXT, verified INTEGER)")
            rows = []
            for name in ("pending.pdf", "new.pdf"):
                (staging / name).write_bytes(b"%PDF-test")
                rows.append(dict(staged_path=name, filename=name, bytes=9))
            md5 = hashlib.md5(b"%PDF-test").hexdigest()
            db.execute("INSERT INTO files VALUES ('pending.pdf', 1, ?, 0)", (md5,)); db.commit()
            pending = dict(id=1, name="pending.pdf", status="ic_checking", supplied_md5=md5)
            client = Mock()
            replies = iter([
                {"location": "file-info"},
                {"id": 2, "status": "created", "upload_url": "upload"},
                {"parts": [{"partNo": 1, "status": "PENDING", "startOffset": 0, "endOffset": 8}]},
                None, None,
            ])
            client.request.side_effect = lambda *a, **kw: next(replies)
            self.assertFalse(transfer_file(client, staging, "record", rows[0], pending))
            client.request.assert_not_called()
            self.assertTrue(transfer_file(client, staging, "record", rows[1], None))
            self.assertEqual(client.request.call_count, 5)
            self.assertEqual(db.execute("SELECT COUNT(*) FROM files WHERE verified=0").fetchone()[0], 2)
            client.request.side_effect = None
            client.request.return_value = {"files": [pending, dict(id=2, name="new.pdf", status="available", computed_md5=md5)]}
            waiting = verify_record(client, db, "record", rows)
            self.assertEqual([r["id"] for r in waiting], [1])
            self.assertEqual(db.execute("SELECT verified FROM files WHERE id=2").fetchone()[0], 1)
            client.request.return_value["files"][1]["computed_md5"] = "wrong"
            with self.assertRaisesRegex(RuntimeError, "Checksum mismatch"):
                verify_record(client, db, "record", rows)
            db.close()

    def test_workers_overlap(self):
        barrier = threading.Barrier(4)
        def task(row):
            barrier.wait(timeout=5)
            return row
        self.assertEqual(sorted(parallel_transfers(task, range(4), 4)), list(range(4)))

    def test_shared_rate_limit_across_workers(self):
        starts = []
        lock = threading.Lock()
        def request(*args, **kwargs):
            with lock:
                starts.append(time.monotonic())
            return Mock(status_code=200, content=b"{}", json=lambda: {})
        session = Mock(request=request)
        with patch("upload_figshare_plots.requests.Session", return_value=session):
            client = Client("test-token")
            list(parallel_transfers(lambda _: client.request("GET", "account"), range(3), 3))
        self.assertTrue(all(b-a >= 1.0 for a, b in zip(starts, starts[1:])))


if __name__ == "__main__":
    unittest.main()
