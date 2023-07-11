import click
from pathlib import Path


@click.command()
@click.option('--bed', '-b', type=Path, required=True,
              help="input bed file")
@click.option('--toml', '-t', type=Path, required=True,
              help="output toml file")
def main(bed, toml):
    with open(toml, 'w') as out:
        out.write('[[jobs]]\n')
        out.write('name = "change_here"\n')
        out.write('time = [0, 6000]\n')
        out.write('ch = [1, 512]\n')
        out.write('bc = "all"\n')
        out.write('target = [\n')

        with open(bed, 'r') as f:
            for line in f:
                content = line.rstrip().split('\t')
                chrom, start, end = content[0], int(content[1]), int(content[2])
                out.write(f'\t{{region = ["{chrom}", "{start}", "{end}"], action = "stop_receiving"}},\n')

        out.write(f'\t{{region = "multi", action = "unblock"}},\n')
        out.write(f'\t{{region = "unmapped", action = "unblock"}},\n')
        out.write(f'\t{{region = "miss", action = "unblock"}},\n')
        out.write("]\n")


if __name__ == "__main__":
    main()
