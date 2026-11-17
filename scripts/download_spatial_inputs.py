"""Download public Xenium inputs listed in the archived input manifest.

Verify both byte counts and SHA256. Never overwrite the archived manifest.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
from urllib.request import Request,urlopen


def digest(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            value.update(block)
    return value.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args()
    folder=args.root/'review_20261013/mural_PDE5A_validation/spatial_GSE318037/input'
    manifest=json.loads((folder/'download_manifest.json').read_text(encoding='utf-8'))
    folder.mkdir(parents=True,exist_ok=True)
    for row in manifest:
        name=row['file']
        if Path(name).name!=name:
            raise ValueError('Manifest filenames must be simple filenames.')
        target=folder/name
        expected=int(row['actual_bytes'])
        expected_hash=row['sha256']
        if target.exists():
            if target.stat().st_size==expected and digest(target)==expected_hash:
                print('Already verified:',name)
                continue
            raise ValueError(f'Existing input has a different size/hash; preserve it for review: {target}')
        partial=target.with_name(target.name+'.partial')
        request=Request(row['url'],headers={'User-Agent':'Public-ocular-spatial-analysis/1.0'})
        with urlopen(request,timeout=120) as response,partial.open('wb') as output:
            shutil.copyfileobj(response,output)
        if partial.stat().st_size!=expected or digest(partial)!=expected_hash:
            raise ValueError(f'Downloaded size/hash mismatch: {partial}')
        partial.replace(target)
        print('Downloaded and verified:',name)


if __name__=='__main__':
    main()
