"""Extract published PDE5A values without renormalization or imputation.

Input preparation only. Existing output is never overwritten.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args()
    source=args.root/'data/raw/GSE183320_all_normalized_expression_data.csv.gz'
    metadata=args.root/'data/raw/GSE183320_all_metadata.csv.gz'
    dest=args.root/'data/processed/sc_choroid_focus_matrix.csv'
    if dest.exists():
        raise FileExistsError(f'Output exists; preserve or explicitly move it before extraction: {dest}')
    if not source.is_file() or not metadata.is_file():
        raise FileNotFoundError('Download both GSE183320 inputs as described in INPUTS.md.')
    found=[]
    for chunk in pd.read_csv(source,chunksize=128):
        gene=chunk.iloc[:,0].astype(str).str.strip()
        for _,row in chunk.loc[gene.eq('PDE5A')].iterrows():
            found.append(row.iloc[1:])
    if len(found)!=1:
        raise ValueError(f'Expected exactly one PDE5A gene row, found {len(found)}.')
    row=found[0]
    values=pd.to_numeric(row,errors='raise').to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values<0).any():
        raise ValueError('Published expression must be finite and nonnegative.')
    result=pd.DataFrame({'cell_barcode':[str(x).strip().replace('-','.') for x in row.index],
                         'PDE5A':values})
    meta=pd.read_csv(metadata,usecols=['barcode'])
    expected=set(meta.barcode.astype(str).str.replace('-','.',regex=False))
    if not result.cell_barcode.is_unique or set(result.cell_barcode)!=expected:
        raise ValueError('Expression cell identifiers do not match the full metadata exactly.')
    dest.parent.mkdir(parents=True,exist_ok=True)
    result.to_csv(dest,index=False)
    print(f'Extracted published PDE5A values for {len(result)} cells to {dest}')


if __name__=='__main__':
    main()
