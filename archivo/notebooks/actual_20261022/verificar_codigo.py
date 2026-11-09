"""Read-only publication integrity, schema and syntax checks; no notebook execution."""
from pathlib import Path
import hashlib
import json
import nbformat

ROOT=Path(__file__).resolve().parent

def verificar():
    manifest=json.loads((ROOT/'manifest_publicacion_20261022.json').read_text(encoding='utf8'))
    for name,digest in manifest['sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    notebooks=list((ROOT/'Notebooks_principales_20261019').glob('*.ipynb'))
    for path in notebooks:
        notebook=nbformat.read(path,as_version=4)
        nbformat.validate(notebook)
        for i,cell in enumerate(notebook.cells):
            if cell.cell_type=='code':
                assert cell.execution_count is None and not cell.outputs,(path.name,i)
                compile(cell.source,f'{path.name}:cell{i}','exec')
    modules=[p for p in ROOT.rglob('*.py') if 'Productos' not in p.parts]
    for path in modules:
        compile(path.read_text(encoding='utf8'),str(path),'exec')
    return {'notebooks':len(notebooks),'modules':len(modules),'integrity':True,'notebooks_executed':False}

if __name__=='__main__':
    print(verificar())
