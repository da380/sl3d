#!/usr/bin/env python3
"""Plot an earth model file.

    ./plot_model.py prem.simp.vm1

The file holds the radius in m, density in kg/m^3, bulk and shear moduli in
GPa, and the base ten logarithm of the viscosity in Pa s, with a zero
viscosity marking a region that is not viscous.
"""

import argparse
import os
import sys

import numpy as np


def read_model(path):
    with open(path) as f:
        name = f.readline().strip()
        nknot, ntau = (int(v) for v in f.readline().split()[:2])
        rows = np.array([[float(v) for v in f.readline().split()[:5]]
                         for _ in range(nknot)])
    return name, ntau, rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file', help='earth model file')
    parser.add_argument('-o', '--output', help='write a png here instead of '
                                               'opening a window')
    args = parser.parse_args()

    import matplotlib
    if args.output:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    name, ntau, rows = read_model(args.file)
    r = rows[:, 0] / 1e3                              # radius in km
    rho, kappa, mu, log_eta = rows[:, 1], rows[:, 2], rows[:, 3], rows[:, 4]

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 5), sharey=True)

    left.plot(rho / 1e3, r, label=r'$\rho$ / 10$^3$ kg m$^{-3}$')
    left.plot(kappa / 100.0, r, label=r'$\kappa$ / 100 GPa')
    left.plot(mu / 100.0, r, label=r'$\mu$ / 100 GPa')
    left.set_xlabel('density and elastic moduli')
    left.set_ylabel('radius / km')
    left.legend(frameon=False)
    left.grid(alpha=0.3)

    # a zero entry marks the elastic lithosphere or the fluid outer core
    viscous = log_eta > 0.0
    right.plot(np.where(viscous, log_eta, np.nan), r, color='C3')
    right.set_xlabel(r'$\log_{10}(\eta$ / Pa s$)$')
    right.grid(alpha=0.3)
    if viscous.any():
        lo, hi = log_eta[viscous].min(), log_eta[viscous].max()
        right.set_xlim(lo - 0.6, hi + 0.6)

    fig.suptitle(f'{os.path.basename(args.file)}'
                 + (f'  ({name}, {ntau} Maxwell element'
                    + ('s' if ntau != 1 else '') + ')' if name else ''))
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f'wrote {args.output}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
