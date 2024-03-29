"""Calculate the covariance of log(rho_t/pi) with log(K(rho_t)(K(pi)))


steps to calculate:
  1. estimate rho_t on grid via one of the methods:
     - histogramming
     - KDE
     - ASH
  2. get pi on same points
  3. Perform convolution with kernel to get K(rho_t) and K(pi)
  4. calculate log of both needed terms for all grid points:
     a = log(rho_t/pi)
     b = log(K(rho_t)/(K(pi)))
  5. put together:
     Cov = \int a*b rho_t dx - (\int a rho_t dx) * (\int b rho_t dx)
"""

import argparse
import numpy as np

from bdld import grid
from bdld.histogram import Histogram
from bdld.potential import polynomial
from bdld.actions.birth_death import dens_kernel_convolution
from helpers.plumed_header import PlumedHeader


def parse_cliargs():
    """Define and get cli arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "trajectories", nargs="+", help="Path to the trajectory files to analyze"
    )
    parser.add_argument(
        "-d",
        "--dim",
        type=int,
        help="Number of dimensions of trajectory, defaults to 1",
        default=1,
    )
    parser.add_argument(
        "--trajectory-stride",
        dest="trajectory_stride",
        type=int,
        help="Stride between trajectory points used as snapshots, defaults to 1",
        default=1,
    )
    parser.add_argument(
        "-kT",
        "--thermal-energy",
        dest="kt",
        type=float,
        help="Energy (in units of kT) of the FES file",
        required=True,
    )
    parser.add_argument(
        "-bw",
        "--birth-death-bandwidth",
        dest="bd_bw",
        nargs="+",
        type=float,
        help="Bandwidth of the birth-death kernel used for the simulation",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--number-of-bins",
        dest="n_bins",
        nargs="+",
        type=int,
        help="Number of bins (for histogramming). "
        "Equal to number of density points being evaluated in KDE case. "
        "Can be omitted for KDE and ASH",
    )
    parser.add_argument(
        "-kde-bw",
        "--kde-bandwidth",
        dest="kde_bw",
        type=float,
        help="KDE bandwidth when using KDE for density estimate",
    )
    parser.add_argument(
        "-m",
        "--method",
        help="Method for density estimation. 'histogram', 'kde' or 'ash'",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Name of the output file, defaults to 'cov_rho_pi'",
        default="cov_rho_pi",
    )
    args = parser.parse_args()

    args.bd_bw = match_dimensions(args.bd_bw, args.dim)
    args.bd_bw = np.array(
        args.bd_bw
    )  # cast to numpy array because some functions expect it

    if args.method == "histogram":
        if not args.n_bins:
            raise ValueError(
                "No number of bins for histogramming specified (--number-of-bins)"
            )
        else:
            args.n_bins = match_dimensions(args.n_bins, args.dim)

    if args.method == "kde":
        if not args.dim == 1:
            raise ValueError("KDE currently only works in 1d")
        if not args.kde_bw:
            raise ValueError("No bandwidth for KDE specified (--kde-bandwidth)")
        if not args.n_bins:  # default for okay resolution
            args.n_bins = [int(2 * 5 / args.kde_bw)]

    return args


def match_dimensions(value_list, dim):
    """Cast single list element to given dimension or check if it matches if more than one"""
    if len(value_list) == 1:
        if dim != 1:  # copy given value to all dimensions
            return value_list * dim
    else:
        if len(value_list) != dim:
            e = "Given number of arguments does not match specified dimensions:\n"
            e += f"{value_list} should have {dim} entries"
            raise ValueError(e)
    return value_list


def main():
    args = parse_cliargs()

    particle_dists, times = read_trajectories(
        args.trajectories, args.dim, args.trajectory_stride
    )

    coeffs = [0, 0.2, -4, 0, 1]
    ranges = [(-2.5, 2.5)]
    pot = polynomial.PolynomialPotential(coeffs, ranges)  # for reference

    # main loop: over the different snapshots at the different times
    values = np.empty((len(particle_dists), 3))  # time, cov, kl_div
    values[:, 0] = times
    for t, dist in enumerate(particle_dists):  # iterate over time
        values[t, 1:] = calc_cov(
            dist,
            pot,
            args.method,
            args.bd_bw,
            args.kt,
            args.n_bins,
            args.kde_bw,
        )

    # set up header for file
    fields = ["time", "cov", "kl_div"]
    constants = {"method": f"{args.method}", "n_bins": f"{args.n_bins}"}
    if args.method == "kde":
        constants["kde_bw"] = f"{args.kde_bw}"
    header = PlumedHeader(fields, constants)

    np.savetxt(
        args.outfile, values, fmt="%14.9f", header=str(header), comments=""
    )


def calc_cov(dist, pot, density_method, bd_bw, kt, n_bins=None, kde_bw_method=None):
    """Calculate the covariance value from a snapshot

    Also returns the "KL divergence" between rho and pi for comparison

    :param dist: numpy array holding the particle distribution
    :param pot: potential used to calculate pi
    :param density_method: method to evaluate density. One of "histogram", "kde", "ash".
    :param bd_bw: bandwidth of the birth-death algorithm during the simulation
    :param kt:  thermal energy
    :param n_bins: List with number of bins for histogram per dimension
                   For KDE this specifies the number of points where the densities
                   are evaluated
    :param kde_bw_method: if set this will use Kernel Density Estimation
                          instead of histogramming. The value is passed to
                          scipy.stats.gaussian_kde (see syntax there)

    :return cov: Full covariance value
    :return kl_div: KL divergence term


    """
    rho, pi = calc_density(dist, pot, density_method, n_bins, kt, kde_bw_method)

    rho_conv = dens_kernel_convolution(rho, bd_bw, "same")
    pi_conv = dens_kernel_convolution(pi, bd_bw, "same")

    a = np.log(rho / pi)
    b = np.log(rho_conv / pi_conv)

    kl_div = integrate(a, rho.data)
    cov = integrate(a * b, rho.data) - kl_div * integrate(b, rho.data)
    return cov, kl_div


def integrate(grid, weights):
    """Integrate grid by simply summing over the rectangles at each point

    The weights are "stronger" than the grid values, i.e. a weight of 0 will
    bring a grid value of "inf" to 0.

    :param grid: grid instance to integrate
    :param weights: weights of the grid points
    :return area: integral value
    """
    area = 1
    for i in range(grid.n_dim):
        area *= (grid.ranges[i][1] - grid.ranges[i][0]) / grid.n_points[i]
    weighted_values = np.where(weights == 0, 0, grid.data * weights)
    return np.sum(weighted_values) * area


def read_trajectories(filenames, dim, stride=1):
    """Read all trajectories from files and return the distributions per time

    :param filenames: list with paths to trajectories
    :param dim: dimensionality of the trajectory data
    :param stride: use only every nth point of the trajectories

    :return particle_dists: numpy array where each element is one distribution snapshot
    :return times: the corresponding times (read from the first file)
    """
    # read in first file to get number of points
    first_dist = np.genfromtxt(filenames[0])[::stride]
    times = first_dist[:, 0]  # assume the simulation time is stored in first column

    # (num_trajectories, number of timesteps, dimensionality of data)
    particle_dists = np.empty((len(filenames), len(times), dim))
    for i, f in enumerate(filenames):
        tmp_dist = np.genfromtxt(f)[::stride]
        if not np.allclose(tmp_dist[:, 0], times):
            raise ValueError(
                f"Simulation times of file {f} do not match the first file"
            )
        particle_dists[i] = tmp_dist[:, 1 : dim + 1]

    # now we have an array with the first index cycling over the trajectories, but we want it to be the time
    particle_dists = np.swapaxes(particle_dists, 0, 1)

    return particle_dists, times


def calc_density(pos, pot, method, n_bins, kt, kde_bw_method=None):
    """Calculate the density (rho) and the corresponding eq. density pi

    this uses either histogramming or kde if

    :param pos: List or numpy array of particle positions
    :param pot: potential used to calculate pi
    :param method: one of "histogram", "kde", "ash"
    :param kt:  thermal energy
    :param n_bins: List with number of bins for histogram per dimension
                   For KDE this specifies the number of points where the densities
                   are evaluated
    :param kde_bw_method: Bandwidth for kde or method to calculate it.
                          The value is passed to scipy.stats.gaussian_kde.

    :return rho: Grid holding the estimated density
    :return pi: Grid with the same points holding the equilibrium density
    """
    if method in ["histogram", "kde"]:  # works similar
        histo = Histogram(n_bins, pot.ranges)
        # rho is the estimated density from the particles
        rho = histo.copy_empty()  # copy points from histo to be consistent between methods
        if method == "histogram":
            histo.add(pos)
            rho.data = histo.data / (np.sum(histo.data) * np.prod(rho.stepsizes))
        elif method == "kde":
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(pos.T, kde_bw_method)
            rho.set_from_func(kde.evaluate)
    elif method == "ash":  # can't use the same grid with the ash implementation
        from ash.ash import ash
        ash_obj = ash(pos.flatten(), n_bins[0])
        rho = grid.from_npoints([(ash_obj.MIN, ash_obj.MAX)], len(ash_obj.ash_mesh))
        rho.data = ash_obj.ash_den
    else:
        raise ValueError(f"Specified method {method} not implemented")

    pot_grid = rho.copy_empty()
    pot_grid.set_from_func(pot.energy)

    # pi is the equilibrium distribution from the potential
    pi = np.exp(-kt * pot_grid)
    pi.data /= np.sum(pi.data) * np.prod(pi.stepsizes)

    return (rho, pi)


if __name__ == "__main__":
    main()
