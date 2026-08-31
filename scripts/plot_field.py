#!/usr/bin/env python3
"""Plot a global field written by one of the forward programs.

    ./plot_field.py sl.64.1d.dacheck.21.0
    ./plot_field.py sl.64.1d.dacheck.21.0 --column 3 -o sea_level.png

Columns are numbered as they are in the file, so --column 4 is the fourth
column, as it would be in gnuplot.  Both forward programs write the sea level
and then its change since the start of the calculation as their last two
columns, which is what --column selects by default, so the same commands work
on gia1d output, which has four columns, and on sl3d output, which also
carries the output index and the projected coordinates.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sl3dplot import read_grid, robinson, draw_coastlines


def find_coastlines(given):
    """Look for the coastline file next to the data, or in the repository."""
    if given:
        return given
    here = os.path.dirname(os.path.abspath(__file__))
    for path in ('coast_lines_110_rob.dat',
                 os.path.join(here, '..', 'data', 'coast_lines_110_rob.dat')):
        if os.path.exists(path):
            return path
    return None


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help='output file from one of the programs')
    parser.add_argument('-c', '--column', type=int,
                        help='column of the file to plot; by default the last, '
                             'which is the change in sea level')
    parser.add_argument('-o', '--output', help='write a png here instead of '
                                               'opening a window')
    parser.add_argument('--cmap', default='RdBu_r', help='matplotlib colour map')
    parser.add_argument('--vmax', type=float,
                        help='saturate the colour scale at this value; by '
                             'default the 99th percentile of the field')
    parser.add_argument('--label',
                        help='label for the colour bar; by default this '
                             'follows the column of the sea level output')
    parser.add_argument('--title', help='title for the plot')
    parser.add_argument('--coastlines', help='projected coastline file')
    args = parser.parse_args()

    import matplotlib
    if args.output:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    lon, lat, columns = read_grid(args.file)
    ncol = columns.shape[0]
    column = args.column if args.column else ncol
    if not 1 <= column <= ncol:
        parser.error(f'{args.file} has columns 1 to {ncol}')
    field = columns[column - 1]

    # close the map by repeating the first meridian at the far edge
    lon = np.append(lon, lon[0] + 360.0)
    field = np.concatenate([field, field[:, :1]], axis=1)

    lon2, lat2 = np.meshgrid(lon, lat)
    x, y = robinson(lat2, lon2)

    vmax = args.vmax if args.vmax else np.percentile(np.abs(field), 99)

    # both forward programs write the sea level and then its change last
    label = args.label or {ncol: 'sea level change / m',
                           ncol - 1: 'sea level / m'}.get(column, 'value')

    fig, ax = plt.subplots(figsize=(11, 6))
    mesh = ax.pcolormesh(x, y, field, cmap=args.cmap, shading='gouraud',
                         vmin=-vmax, vmax=vmax)

    coast = find_coastlines(args.coastlines)
    if coast:
        draw_coastlines(ax, coast)
    else:
        print('no coastline file found; plotting without it', file=sys.stderr)

    # the outline of the projection
    edge_lat = np.linspace(-90, 90, 181)
    for edge_lon in (-180.0, 180.0):
        ex, ey = robinson(edge_lat, np.full_like(edge_lat, edge_lon))
        ax.plot(ex, ey, color='black', linewidth=0.8)
    for edge in (-90.0, 90.0):
        ex, ey = robinson(np.full(2, edge), np.array([-180.0, 180.0]))
        ax.plot(ex, ey, color='black', linewidth=0.8)

    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.set_title(args.title if args.title
                 else f'{os.path.basename(args.file)}, column {column}')
    fig.colorbar(mesh, ax=ax, orientation='horizontal', shrink=0.6, pad=0.02,
                 label=label)
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches='tight')
        print(f'wrote {args.output}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
