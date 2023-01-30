#!/usr/bin/env python3
import click
import pandas as pd


@click.command()
@click.option('--fl', type=str, required=True,
              help='quant.sf for full-length reads.')
@click.option('--nonfl', type=str, required=True,
              help='quant.sf for non-fl reads.')
@click.option('--output', type=str, required=True,
              help='output quantification result.')
def main(fl, nonfl, output):
    col = "NumReads"
    fl_sf = pd.read_csv(fl, sep="\t", index_col=0)
    nonfl_sf = pd.read_csv(nonfl, sep="\t", index_col=0)

    merged_reads = fl_sf[col] + nonfl_sf[col]
    merged_sf = pd.DataFrame({
        "Length": fl_sf['Length'],
        "EffectiveLength": fl_sf['EffectiveLength'],
        "CPM": 1000 * 1000 * merged_reads / merged_reads.sum(),
        "NumReads": merged_reads}
    )
    merged_sf.to_csv(output, sep="\t", index=True, index_label="Name")


if __name__ == '__main__':
    main()