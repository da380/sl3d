"""Helpers shared by the plotting scripts.

Reading the output grids the codes write, the Robinson projection they use,
and the coastlines that go with it.  Only numpy is needed here; matplotlib is
used by the scripts themselves.
"""

import numpy as np


# ---------------------------------------------------------------- output files


def read_grid(path):
    """Read one of the codes' output files.

    The forward programs write a global field as one row per grid point, in
    blocks of constant longitude separated by blank lines, which is the layout
    gnuplot's pm3d wants.  What the columns hold differs between the two codes:
    gia1d writes longitude, latitude, sea level and its change since the start
    of the calculation, while sl3d writes the output index and the projected
    coordinates as well.  The blocks give away the shape of the grid, and the
    latitude is the first column that varies inside a block, so both layouts
    can be read without being told which is which.

    Returns (lon, lat, columns), where columns has shape (ncol, nlat, nlon) and
    is indexed by column of the file, so columns[0] is the first column.
    """
    blocks, current = [], []
    with open(path) as f:
        for line in f:
            if line.strip():
                current.append([float(v) for v in line.split()])
            elif current:
                blocks.append(current)
                current = []
    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError(f'{path}: no data')

    width = len(blocks[0])
    blocks = [b for b in blocks if len(b) == width]
    grid = np.array(blocks)                       # (nlon, nlat, ncol)
    if grid.ndim != 3 or grid.shape[2] < 3:
        raise ValueError(f'{path}: expected at least three columns')

    # latitude varies down a block; longitude is the column before it
    varies = [c for c in range(grid.shape[2])
              if not np.allclose(grid[0, :, c], grid[0, 0, c])]
    if not varies or varies[0] == 0:
        raise ValueError(f'{path}: cannot tell which columns are the grid')
    ilat = varies[0]
    ilon = ilat - 1

    lon = grid[:, 0, ilon]
    lat = grid[0, :, ilat]
    return lon, lat, np.transpose(grid, (2, 1, 0))


def sample(lon, lat, field, site_lat, site_lon):
    """Bilinear sample of a field at one point, wrapping in longitude."""
    j = np.interp(site_lon, lon, np.arange(lon.size), period=360.0)
    i = np.interp(site_lat, lat, np.arange(lat.size))
    i0, j0 = int(np.floor(i)), int(np.floor(j))
    i1, j1 = min(i0 + 1, lat.size - 1), (j0 + 1) % lon.size
    di, dj = i - i0, j - j0
    return ((1 - di) * (1 - dj) * field[i0, j0] + (1 - di) * dj * field[i0, j1]
            + di * (1 - dj) * field[i1, j0] + di * dj * field[i1, j1])


# ------------------------------------------------------------ map projection


# The Robinson table, at five degree intervals from the equator to the pole,
# exactly as it appears in module_gia.f90.
_PLEN = np.array([1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600,
                  0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597,
                  0.7186, 0.6732, 0.6213, 0.5722, 0.5322])
_PDFE = np.array([0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720,
                  0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903,
                  0.8435, 0.8936, 0.9394, 0.9761, 1.0000])
_LATS = np.concatenate([-5.0 * np.arange(18, 0, -1), 5.0 * np.arange(0, 19)])
_TABLE = {
    'plen': np.concatenate([_PLEN[:0:-1], _PLEN]),
    'pdfe': np.concatenate([-_PDFE[:0:-1], _PDFE]),
}


def _natural_spline(x, y):
    """Second derivatives of the natural cubic spline through (x, y).

    The Fortran interpolates its Robinson table with a natural cubic spline,
    and doing the same here keeps a plotted field registered with the
    coastlines in data/coast_lines_110_rob.dat, which were projected by it.
    """
    n = x.size
    h = np.diff(x)
    a = np.zeros((n, n))
    r = np.zeros(n)
    a[0, 0] = a[-1, -1] = 1.0                       # natural: y'' = 0 at both ends
    for i in range(1, n - 1):
        a[i, i - 1] = h[i - 1] / 6.0
        a[i, i] = (h[i - 1] + h[i]) / 3.0
        a[i, i + 1] = h[i] / 6.0
        r[i] = (y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]
    return np.linalg.solve(a, r)


def _splint(x, y, y2, xi):
    """Evaluate the spline defined by _natural_spline."""
    xi = np.asarray(xi, dtype=float)
    k = np.clip(np.searchsorted(x, xi) - 1, 0, x.size - 2)
    h = x[k + 1] - x[k]
    a = (x[k + 1] - xi) / h
    b = (xi - x[k]) / h
    return (a * y[k] + b * y[k + 1]
            + ((a ** 3 - a) * y2[k] + (b ** 3 - b) * y2[k + 1]) * h ** 2 / 6.0)


_COEF = {name: _natural_spline(_LATS, tab) for name, tab in _TABLE.items()}


def robinson(lat, lon):
    """Project latitude and longitude in degrees, as get_xy does in the Fortran.

    The map spans roughly -0.5 to 0.5 in x and -0.507 to 0.507 in y.
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    plen = _splint(_LATS, _TABLE['plen'], _COEF['plen'], lat)
    pdfe = _splint(_LATS, _TABLE['pdfe'], _COEF['pdfe'], lat)
    return lon * plen / 360.0, pdfe * 0.5072


def coastlines(path):
    """Read the projected coastlines, returning a list of line segments."""
    segments, current = [], []
    with open(path) as f:
        for line in f:
            if line.strip():
                x, y = line.split()[:2]
                current.append((float(x), float(y)))
            elif current:
                segments.append(np.array(current))
                current = []
    if current:
        segments.append(np.array(current))
    return segments


def draw_coastlines(ax, path, **kwargs):
    """Draw the coastlines onto a set of axes."""
    style = dict(color='black', linewidth=0.6)
    style.update(kwargs)
    for seg in coastlines(path):
        ax.plot(seg[:, 0], seg[:, 1], **style)
