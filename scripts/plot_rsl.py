#!/usr/bin/env python3
"""Plot a relative sea level curve at one site from a run's output files.

    ./plot_rsl.py --lat 57.7 --lon 11.9 sl.64.1d.dacheck.*

The forward programs write one file per output time, with the age in thousands
of years at the end of the name, so the whole history of a site can be read
off a run's output.  Several sites may be given by repeating --site.
"""

import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sl3dplot import read_grid, sample


def age_of(path):
    """The age in ka from the end of a file name, as the programs write it."""
    m = re.search(r'(\d+\.\d+)$', os.path.basename(path))
    return float(m.group(1)) if m else None


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('files', nargs='+', help='output files of one run')
    parser.add_argument('--site', action='append', nargs=3,
                        metavar=('LAT', 'LON', 'NAME'),
                        help='a site to plot; may be repeated')
    parser.add_argument('--lat', type=float, help='latitude of a single site')
    parser.add_argument('--lon', type=float, help='longitude of a single site')
    parser.add_argument('-c', '--column', type=int,
                        help='column of the files to plot; by default the '
                             'last, which is the change in sea level')
    parser.add_argument('-o', '--output', help='write a png here instead of '
                                               'opening a window')
    args = parser.parse_args()

    sites = [(float(a), float(b), c) for a, b, c in (args.site or [])]
    if args.lat is not None and args.lon is not None:
        sites.append((args.lat, args.lon, f'{args.lat:.1f}, {args.lon:.1f}'))
    if not sites:
        parser.error('give a site with --lat and --lon, or with --site')

    dated = [(age_of(f), f) for f in args.files]
    if any(age is None for age, _ in dated):
        parser.error('every file name must end in an age, as in '
                     'sl.64.1d.dacheck.20.5')
    dated.sort(reverse=True)          # oldest first

    import matplotlib
    if args.output:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ages = np.array([age for age, _ in dated])
    curves = {name: [] for _, _, name in sites}
    column = None
    for _, path in dated:
        lon, lat, columns = read_grid(path)
        ncol = columns.shape[0]
        column = args.column if args.column else ncol
        if not 1 <= column <= ncol:
            parser.error(f'{path} has columns 1 to {ncol}')
        field = columns[column - 1]
        for site_lat, site_lon, name in sites:
            curves[name].append(sample(lon, lat, field, site_lat, site_lon))

    fig, ax = plt.subplots(figsize=(8, 5))
    for name, values in curves.items():
        ax.plot(ages, values, marker='o', markersize=3, label=name)
    ax.invert_xaxis()                 # age runs backwards to the present
    ax.axhline(0.0, color='grey', linewidth=0.8, zorder=0)
    ax.set_xlabel('age / ka')
    ax.set_ylabel({ncol: 'sea level change / m',
                   ncol - 1: 'sea level / m'}.get(column, 'value'))
    ax.set_title('relative sea level')
    ax.legend(frameon=False)
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f'wrote {args.output}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
