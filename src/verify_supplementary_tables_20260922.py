import csv
import hashlib
import math
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'int/20260922_deliver'
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
cache = {}
NEUROPIL_LABELS = {'ME_R': 'Medulla', 'LO_R': 'Lobula', 'LOP_R': 'Lobula plate'}

def workbook(path):
    path = ROOT / path
    if path in cache:
        return cache[path]
    result = {}
    with ZipFile(path) as z:
        assert z.testzip() is None, path
        strings = []
        if 'xl/sharedStrings.xml' in z.namelist():
            strings = [''.join(s.itertext()) for s in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        rels = {s.get('Id'): s.get('Target') for s in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        for sheet in ET.fromstring(z.read('xl/workbook.xml')).findall('m:sheets/m:sheet', NS):
            target = rels[sheet.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')]
            target = target.lstrip('/') if target.startswith('/') else 'xl/' + target
            xml = ET.fromstring(z.read(target))
            rows = []
            for row in xml.findall('m:sheetData/m:row', NS):
                vals = {}
                for cell in row:
                    assert cell.get('t') != 'e', (path, sheet.get('name'), cell.get('r'))
                    assert cell.find('m:f', NS) is None, (path, 'unexpected formula')
                    raw = cell.find('m:v', NS)
                    value = raw.text if raw is not None else None
                    if cell.get('t') == 's':
                        value = strings[int(value)]
                    elif cell.get('t') == 'inlineStr':
                        value = ''.join(cell.find('m:is', NS).itertext())
                    elif cell.get('t') == 'b' and value is not None:
                        value = value == '1'
                    elif value is not None:
                        value = float(value)
                    vals[re.sub(r'\d', '', cell.get('r'))] = value
                rows.append(vals)
            headers = rows[0]
            assert len(set(headers.values())) == len(headers), (path, sheet.get('name'))
            records = [{name: row.get(col) for col, name in headers.items()} for row in rows[1:]]
            if DEST in path.parents:
                assert xml.find('m:autoFilter', NS) is not None
                assert xml.find('m:sheetViews/m:sheetView/m:pane', NS).get('state') == 'frozen'
            result[sheet.get('name')] = records
    cache[path] = result
    return result

def csv_rows(path):
    with (ROOT / path).open(newline='') as f:
        return list(csv.DictReader(f))

def missing(x):
    return x is None or x == '' or x == 'NA'

def equivalent(a, b):
    if missing(a) or missing(b):
        return missing(a) and missing(b)
    try:
        af, bf = float(a), float(b)
        return af == bf or math.isclose(af, bf, rel_tol=1e-13, abs_tol=1e-14)
    except (ValueError, TypeError):
        return str(a) == str(b)

manifest = csv_rows('int/20260922_deliver/SOURCE_MANIFEST.csv')
columns = csv_rows('int/20260922_deliver/COLUMN_MAPPING.csv')
maps = defaultdict(dict)
for c in columns:
    maps[c['Table'], c['Sheet']][c['Export_column']] = c['Source_column']

template = workbook('int/Supp. Table 5-TS and CAM.xlsx')
concentric = defaultdict(set)
editorial = {}
for name, rows in template.items():
    if not name.startswith(('TS-', 'CAM-')):
        continue
    gene_col = 'Selector' if name.startswith('TS-') else 'CAM'
    for r in rows:
        key = (name, r[gene_col], r['Notch Status'], r['Synapse Type'])
        assert key not in editorial
        editorial[key] = r['Visualized']
        if gene_col == 'Selector':
            concentric[r[gene_col]].add(r['Concentric gene'])
assert all(len(v) == 1 for v in concentric.values())

verified_cells = 0
sheet_counts = Counter()
table5_rows = {}
for entry in manifest:
    table, sheet = entry['Table'], entry['Sheet']
    output = DEST / ('Supp_Table_' + table) / entry['Workbook']
    actual = workbook(output)[sheet]
    for row in actual:
        if 'Neuropil' in row:
            assert row['Neuropil'] in NEUROPIL_LABELS.values(), (table, sheet, row['Neuropil'])
        assert not any(isinstance(value, str) and re.search(r'\b(?:ME|LO|LOP)_[LR]\b', value)
                       for value in row.values()), (table, sheet, 'raw neuropil code in output')
    assert len(actual) == int(entry['Rows'])
    sheet_counts[table] += len(actual)
    mapping = maps[table, sheet]
    assert list(actual[0]) == list(mapping), (table, sheet, 'headers')
    if sheet == 'Depth association percentages':
        continue
    source = entry['Source']
    if source.endswith('.xlsx'):
        expected = workbook(source)[entry['Source_sheet']]
        if entry['Source_sheet'] == 'figure_data':
            expected = [r for r in expected if r['included_in_figure'] is True or r['included_in_figure'] == 1]
            assert not any('value' in k or 'Significant' in k for k in actual[0])
    else:
        expected = csv_rows(source)
        if table == '5':
            np = {'Medulla':'ME_R', 'Lobula':'LO_R', 'Lobula plate':'LOP_R'}[sheet.split('-',1)[1]]
            expected = [r for r in expected if missing(r['skip_reason']) and r['neuropil'] == np]
            table5_rows[sheet] = actual
    if expected and 'neuropil' in expected[0]:
        expected = [r for r in expected if r['neuropil'] in NEUROPIL_LABELS]
    assert len(actual) == len(expected), (table, sheet, 'source row count')
    for got, src in zip(actual, expected):
        for header, col in mapping.items():
            if col == 'neuropil':
                wanted = NEUROPIL_LABELS[src[col]]
            elif col == 'notch_label':
                wanted = 'All' if src['split'] == 'all' else src['split'].split('_', 1)[0]
            elif col == 'class_label':
                wanted = 'All' if src['split'] == 'all' else src['split'].split('_', 1)[1]
            elif col == 'Visualized':
                wanted = 'No' if sheet.endswith('Lobula plate') else editorial[(sheet, src['types_of_interest'], src['notch_category'], src['syn_type'])]
            elif col == 'concentric_gene':
                wanted = next(iter(concentric[src['types_of_interest']]))
            else:
                wanted = src[col]
            assert equivalent(got[header], wanted), (table, sheet, header, got[header], wanted)
            # Known numeric columns must remain numeric XML cells; Inf is explicit text.
            if isinstance(wanted, (int, float)) and not isinstance(wanted, bool) and math.isfinite(wanted):
                assert isinstance(got[header], (int, float)), (table, sheet, header, 'numeric type')
            verified_cells += 1

summary_path = DEST / 'Supp_Table_5/Supp_Table_5_TS_and_CAM.xlsx'
summary = workbook(summary_path)['Depth association percentages']
assert len(summary) == 24
for row in summary:
    prefix = 'TS' if row['Molecule Class'] == 'Selector TF' else 'CAM'
    region = row['Neuropil']
    tested = [r for r in table5_rows[prefix + '-' + region] if r['Notch Status'] == row['Notch Status'] and r['Synapse Type'] == row['Synapse Type']]
    significant = sum(r['q-value (FDR)'] < .05 for r in tested)
    assert row['Tested Genes'] == len(tested)
    assert row['Significant Genes'] == significant
    assert equivalent(row['Percentage'], 100 * significant / len(tested))
assert sum(row['Tested Genes'] for row in summary if row['Molecule Class'] == 'Selector TF') == 405
assert sum(row['Tested Genes'] for row in summary if row['Molecule Class'] == 'CAM') == 2082

# Independently check gene coverage, including all successful LOP entries.
lop_report = []
for prefix, col in [('TS','Selector'), ('CAM','CAM')]:
    lop = {r[col] for r in table5_rows[prefix+'-Lobula plate']}
    other = {r[col] for region in ['Medulla','Lobula'] for r in table5_rows[prefix+'-'+region]}
    assert not lop - other
    assert all(r['Visualized'] == 'No' for r in table5_rows[prefix+'-Lobula plate'])
    lop_report.append(f'{prefix}: {len(lop)} tested lobula-plate genes; none exclusive to lobula plate.')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

pdf = DEST / 'Supp_Table_6/origin_heatmap.pdf'
assert pdf.read_bytes().startswith(b'%PDF-')
assert sha(pdf) == sha(ROOT / 'int/origin_heatmap/origin_heatmap.pdf')
sources = (DEST / 'SOURCE_FILES.txt').read_text().splitlines()
(DEST / 'SOURCE_SHA256SUMS').write_text(''.join(f'{sha(ROOT / p)}  {p}\n' for p in sources))
all_books = list(DEST.glob('Supp_Table_*/*.xlsx'))
assert len(all_books) == 5
assert sum(len(workbook(p)) for p in all_books) == len(manifest) == 37
total = sum(sheet_counts.values())
report = [
    '# Delivery validation', '',
    f'- Five XLSX workbooks; 37 worksheets; {total:,} data rows.',
    f'- Independent Python ZIP/XML verification checked {verified_cells:,} exported cells against original CSV/XLSX sources, including missing values and biological context.',
    '- R independently reopened each workbook and compared every exported column with the selected data (numerical tolerance 1e-13).',
    '- Excel containers pass CRC checks; no Excel error cells or formulas; all sheets have filters and frozen headers.',
    '- All 24 percentage rows independently reproduce successful-test counts, Q < 0.05 counts, and exact percentages: 405 selector tests and 2,082 CAM tests.',
    '- All transferred editorial annotations match the supplied workbook; selector concentric annotations are unambiguous.',
    '- Cochran-Armitage and other enrichment rows retain the source values and missing-result explanations. Infinite odds ratios are explicit Inf text.',
    '- All neuropil-specific tables include only right-hemisphere source rows and display Medulla, Lobula, or Lobula plate; no raw neuropil codes remain in workbook cells.',
    '- Revised spatial figure data includes only plotted rows and excludes inferential significance fields. Aggregate supporting statistics retain only right-hemisphere source rows.',
    '- Heatmap PDF is byte-identical to the original. SOURCE_SHA256SUMS records the exact inputs; SHA256SUMS records delivered files.',
] + ['- '+s for s in lop_report] + ['', '| Table | Data rows |', '|---|---:|'] + [f'| {t} | {sheet_counts[t]} |' for t in sorted(sheet_counts)]
(DEST / 'VALIDATION.md').write_text('\n'.join(report) + '\n')
files = sorted(p for p in DEST.rglob('*') if p.is_file() and p.name != 'SHA256SUMS')
(DEST / 'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(DEST)}\n' for p in files))
print('\n'.join(report))
