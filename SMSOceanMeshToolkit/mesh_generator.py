# stdlib
import logging
import time

# tpl
import numpy as np
import scipy.sparse as spsparse
import matplotlib.pyplot as plt

from .clean import fix_mesh, simp_qual
from .Grid import Grid
from .libs._delaunay_class import DelaunayTriangulation as DT
from .libs._fast_geometry import unique_edges
# local
from .signed_distance_function import SDFDomain as Domain
from .plotting import SimplexCollection

logger = logging.getLogger(__name__)

__all__ = ["generate_mesh"]

DENSITY_CONTROL_FREQUENCY = 30

def _parse_kwargs(kwargs):
    for key in kwargs:
        if key in {
            "nscreen",
            "max_iter",
            "seed",
            "pfix",
            "points",
            "domain",
            "edge_length",
            "bbox",
            "min_edge_length",
            "plot",
            "pseudo_dt",
        }:
            pass
        else:
            raise ValueError(
                "Option %s with parameter %s not recognized " % (key, kwargs[key])
            )


def _check_bbox(bbox):
    assert isinstance(bbox, tuple), "`bbox` must be a tuple"
    assert int(len(bbox) / 2), "`dim` must be 2"

def _plot_mesh(fig, c, bbox, p, t, count): 
    xmin, xmax = bbox[0]
    ymin, ymax = bbox[1]
    # red to green colormap spanning 11 colors 
    if count==0:
        fig.clf()
        ax = fig.gca()
        c = SimplexCollection()
        ax.add_collection(c)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal')
        ax.set_axis_off()
        c.set_simplices((p, t), color_by_quality=True)
        fig.canvas.draw()    
    else:   
        c.set_simplices((p, t), color_by_quality=True)
        fig.canvas.draw()
    # set the title to the count 
    fig.suptitle(f"Iteration {count}")
    return fig, c 

        


def generate_mesh(domain, edge_length, **kwargs):
    """
    Generate a 2D triangular mesh using callbacks to a
    mesh sizing function `edge_length` and signed distance function `domain`

    Parameters
    ----------
    domain: A function object.
        A function that takes a point and returns the signed nearest distance to the domain boundary Ω.
    edge_length: A function object.
        A function that can evalulate a point and return a mesh size.
    \**kwargs:
        See below

    :Keyword Arguments:
        * *bbox* (``tuple``) --
            Bounding box containing domain extents. REQUIRED IF NOT USING :class:`edge_length`
        * *max_iter* (``float``) --
            Maximum number of meshing iterations. (default==50)
        * *seed* (``float`` or ``int``) --
            Pseudo-random seed to initialize meshing points. (default==0)
        * *pfix* (`array-like`) --
            An array of points to constrain in the mesh. (default==None)
        * *min_edge_length* (``float``) --
            The minimum element size in the domain. REQUIRED IF NOT USING :class:`edge_length`
        * *plot* (``int``) --
            The mesh is visualized every `plot` meshing iterations.
        * *pseudo_dt* (``float``) --
            The pseudo time step for the meshing algorithm. (default==0.2)

    Returns
    -------
    points: array-like
        vertex coordinates of mesh
    t: array-like
        mesh connectivity table.

    """
    _DIM = 2
    opts = {
        "max_iter": 50,
        "seed": 0,
        "pfix": None,
        "points": None,
        "min_edge_length": None,
        "plot": 999999,
        "pseudo_dt": 0.2,
    }
    opts.update(kwargs)
    _parse_kwargs(kwargs)

    fd, bbox = _unpack_domain(domain, opts)
    fh, min_edge_length = _unpack_sizing(edge_length, opts)

    _check_bbox(bbox)
    bbox = np.array(bbox).reshape(-1, 2)

    assert min_edge_length > 0, "`min_edge_length` must be > 0"

    assert opts["max_iter"] > 0, "`max_iter` must be > 0"
    max_iter = opts["max_iter"]
    # log the seed used 
    logger.info(f"Using seed {opts['seed']}")
    np.random.seed(opts["seed"])

    L0mult = 1 + 0.4 / 2 ** (_DIM - 1)
    delta_t = opts["pseudo_dt"]
    #geps = 1e-3 * np.amin(min_edge_length)
    Re = 6378.137e3
    geps = 1e-12*min_edge_length/Re;
    deps = np.sqrt(np.finfo(np.double).eps) * np.amin(min_edge_length)

    pfix, nfix = _unpack_pfix(_DIM, opts)

    if opts["points"] is None:
        logging.info("Generating initial points")
        p = _generate_initial_points(
            min_edge_length,
            geps,
            bbox,
            fh,
            fd,
            pfix,
        )
    else:
        logging.info("Using user-provided points")
        p = opts["points"]
        # assert that it's 2d 
        assert p.shape[1] == 2, "Points must be 2D"

    N = p.shape[0]

    assert N > 0, "No vertices to mesh with!"

    logger.info(
        f"Commencing mesh generation with {N} vertices will perform {max_iter} iterations."
    )
    if opts["plot"] > 0:
        fig = plt.figure()
        
    for count in range(max_iter):
        start = time.time()

        # (Re)-triangulation by the Delaunay algorithm
        dt = DT()
        dt.insert(p.ravel().tolist())

        # Get the current topology of the triangulation
        p, t = _get_topology(dt)

        ifix = []
        # Find where pfix went
        if nfix > 0:
            for fix in pfix:
                ind = _closest_node(fix, p)
                ifix.append(ind)
                p[ind] = fix

        # Remove points outside the domain
        t = _remove_triangles_outside(p, t, fd, geps)

        if count % opts["plot"] == 0:
            if count == 0:
                c = None
            fig, c = _plot_mesh(fig, c, bbox, p, t, count)
            # For visualizing the mesh
            plt.pause(0.2)

        # Number of iterations reached, stop.
        if count == (max_iter - 1):
            p, t, _ = fix_mesh(p, t, dim=_DIM, delete_unused=True)
            logger.info("Termination reached...maximum number of iterations.")
            plt.close() # close any plots
            return p, t

        # Compute the forces on the bars
        Ftot, p, N, pold = _compute_forces(p, t, fh, L0mult, nfix, count)
        if np.all(Ftot == 0):
            continue
        #Ftot = _compute_forces(p, t, fh, min_edge_length, L0mult)
        
        # Force = 0 at fixed points
        Ftot[:nfix] = 0

        # Update positions
        p += delta_t * Ftot

        # Bring outside points back to the boundary
        p = _project_points_back(p, fd, deps)

        # Show the user some progress so they know something is happening
        maxdp = delta_t * np.sqrt((Ftot**2).sum(1)).max()

        logger.info(
            f"Iteration #{count + 1}, max movement is {maxdp}, there are {len(p)} vertices and {len(t)}"
        )

        end = time.time()
        logger.info(f"Elapsed wall-clock time {end - start} seconds")


def _unpack_sizing(edge_length, opts):
    if isinstance(edge_length, Grid):
        fh = edge_length.eval
        min_edge_length = np.min(edge_length.values)
    elif callable(edge_length):
        fh = edge_length
        min_edge_length = opts["min_edge_length"]
    else:
        raise ValueError(
            "`edge_length` must either be a function or a `edge_length` object"
        )
    return fh, min_edge_length


def _unpack_domain(domain, opts):
    if isinstance(domain, Domain):
        bbox = domain.bbox
        fd = domain.eval
    elif callable(domain):
        bbox = opts["bbox"]
        fd = domain
    else:
        raise ValueError(
            "`domain` must be a function or a :class:`signed_distance_function object"
        )
    return fd, bbox


def _get_bars(t):
    """Describe each bar by a unique pair of nodes"""
    bars = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    return unique_edges(bars)


# Persson-Strang
def _compute_forces(p, t, fh, L0mult, nfix, count):
    """Compute the forces on each edge based on the sizing function"""
    N = p.shape[0]
    bars = _get_bars(t)
    barvec = p[bars[:, 0]] - p[bars[:, 1]]  # List of bar vectors
    L = np.sqrt((barvec**2).sum(1))  # L = Bar lengths
    L[L == 0] = np.finfo(float).eps
    hbars = fh(p[bars].sum(1) / 2)
    # if hbars is excessively large then abort 
    if np.any(hbars > 1e4):
        raise ValueError("Problem detected in fh...is extrapolation enabled?")
    L0 = hbars * L0mult * (np.nanmedian(L) / np.nanmedian(hbars)) 
    F = L0 - L
    F[F < 0] = 0  # Bar forces (scalars)
    Fvec = (
        F[:, None] / L[:, None].dot(np.ones((1, 2))) * barvec
    )  # Bar forces (x,y components)
    Ftot = _dense(
        bars[:, [0] * 2 + [1] * 2],
        np.repeat([list(range(2)) * 2], len(F), axis=0),
        np.hstack((Fvec, -Fvec)),
        shape=(N, 2),
    )
    pold = p 
    # Density control - remove points that are too close
    if (count % DENSITY_CONTROL_FREQUENCY) == 0 and (L0 > 2*L).any():
        ixdel = np.setdiff1d(bars[L0 > 2*L].reshape(-1), np.arange(nfix))
        p = p[np.setdiff1d(np.arange(N), ixdel)]
        N = p.shape[0]; pold = float('inf')
        logger.info(f"Removed {len(ixdel)} points due to excessive density")
        Ftot *= 0.0 # reset the forces
        
    return Ftot, p, N, pold


# Bossen-Heckbert
#def _compute_forces(p, t, fh, min_edge_length, L0mult):
#   """Compute the forces on each edge based on the sizing function"""
#   N = p.shape[0]
#   bars = _get_bars(t)
#   barvec = p[bars[:, 0]] - p[bars[:, 1]]  # List of bar vectors
#   L = np.sqrt((barvec ** 2).sum(1))  # L = Bar lengths
#   L[L == 0] = np.finfo(float).eps
#   hbars = fh(p[bars].sum(1) / 2)
#   L0 = hbars * L0mult * (np.nanmedian(L) / np.nanmedian(hbars))
#   LN = L / L0
#   F = (1 - LN ** 4) * np.exp(-(LN ** 4)) / LN
#   Fvec = (
#       F[:, None] / LN[:, None].dot(np.ones((1, 2))) * barvec
#   )  # Bar forces (x,y components)
#   Ftot = _dense(
#       bars[:, [0] * 2 + [1] * 2],
#       np.repeat([list(range(2)) * 2], len(F), axis=0),
#       np.hstack((Fvec, -Fvec)),
#       shape=(N, 2),
#   )
#   return Ftot


def _dense(Ix, J, S, shape=None, dtype=None):
    """
    Similar to MATLAB's SPARSE(I, J, S, ...), but instead returning a
    dense array.
    """

    # Advanced usage: allow J and S to be scalars.
    if np.isscalar(J):
        x = J
        J = np.empty(Ix.shape, dtype=int)
        J.fill(x)
    if np.isscalar(S):
        x = S
        S = np.empty(Ix.shape)
        S.fill(x)

    # Turn these into 1-d arrays for processing.
    S = S.flat
    II = Ix.flat
    J = J.flat
    return spsparse.coo_matrix((S, (II, J)), shape, dtype).toarray()


def _remove_triangles_outside(p, t, fd, geps):
    """Remove vertices outside the domain"""
    pmid = p[t].sum(1) / 3  # Compute centroids
    return t[fd(pmid) < -geps]  # Keep interior triangles


def _project_points_back(p, fd, deps):
    """Project points outsidt the domain back within"""
    d = fd(p)
    ix = d > 0  # Find points outside (d>0)
    if ix.any():

        def _deps_vec(i):
            a = [0] * 2
            a[i] = deps
            return a

        try:
            dgrads = [
                (fd(p[ix] + _deps_vec(i)) - d[ix]) / deps for i in range(2)
            ]  # old method
        except ValueError:  # an error is thrown if all points in fd are outside
            # bbox domain, so instead calulate all fd and then
            # take the solely ones outside domain
            dgrads = [(fd(p + _deps_vec(i)) - d) / deps for i in range(2)]
            dgrads = list(np.array(dgrads)[:, ix])
        dgrad2 = sum(dgrad**2 for dgrad in dgrads)
        dgrad2 = np.where(dgrad2 < deps, deps, dgrad2)
        p[ix] -= (d[ix] * np.vstack(dgrads) / dgrad2).T  # Project
    return p

def _generate_initial_points(min_edge_length, geps, bbox, fh, fd, pfix):
    """
    Generate an initial distribution of points within a bounding box based on equilateral triangles.

    Parameters:
    - bbox: Bounding box as a list of tuples [(min_x, max_x), (min_y, max_y)].
    - fh: Function handler to adjust point density.
    - fd: Function handler to define the domain, points where fd < 0 will be kept.
    - min_edge_length: Minimum edge length for equilateral triangles.
    - geps: Geometry epsilon, tolerance for defining inside/outside of the domain.
    - pfix: Array of fixed points that must be included in the output.

    Returns:
    - A numpy array of points within the domain defined by fd, adjusted by fh, including pfix.
    """
    xmin, xmax, ymin, ymax = bbox.flatten()
    # create a polygon from the bbox 
    from shapely.geometry import Polygon
    bbox_tmp = Polygon([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)])
    
    h0 = min_edge_length
    # Create initial distribution in bounding box (equilateral triangles)
    x, y = np.mgrid[xmin:(xmax+h0):h0,
                    ymin:(ymax+h0*np.sqrt(3)/2):h0*np.sqrt(3)/2]
    #x[:, 1::2] += h0/2                               # Shift even rows
    y[:, 1::2] += h0*np.sqrt(3)/4                     # Shift even rows
    p = np.vstack((x.flat, y.flat)).T                # List of node coordin
    
    # Filter points using fh and fd
    p = p[fd(p) < geps]
    r0 = 1/fh(p)**2                                  # Probability to keep point
    p = p[np.random.random(p.shape[0])<r0/r0.max()]  # Rejection method

    # Concatenate fixed points
    if pfix is not None and len(pfix) > 0:
        p = np.vstack([pfix, p])

    return p

def _unpack_pfix(dim, opts):
    """Unpack fixed points"""
    pfix = np.empty((0, dim))
    nfix = 0
    if opts["pfix"] is not None:
        pfix = np.array(opts["pfix"], dtype="d")
        nfix = len(pfix)
        logger.info(f"Constraining {nfix} fixed points..")
    return pfix, nfix


def _get_topology(dt):
    """Get points and entities from :clas:`CGAL:DelaunayTriangulation2/3` object"""
    return dt.get_finite_vertices(), dt.get_finite_cells()


def _closest_node(node, nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist_2 = np.einsum("ij,ij->i", deltas, deltas)
    return np.argmin(dist_2)
