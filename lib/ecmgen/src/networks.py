from .stranddistributions import (
    UniformStrandDistribution,
    DeterministicStrandDistribution,
    StrandDistributionGeneral,
    VonMisesStrandDistribution,
)
from .strandgens import RandomStrandGenerator
from .crosslink_distributors import TipToTailCrosslinkDistributer
from .density_crosslinker import StrandDensityCrosslinkDistributer
from .parameters import (
    DomainParameters,
    RandomStrandGeneratorParameters,
    StrandDensityCrosslinkDistributerParameters,
)
from .density_crosslinker import (
    StrandDensityCrosslinkDistributer,
    StrandDensityCrosslinkDistributerFast,
)
from .crosslink_distributors import DeterministicCrosslinkDistributer
from .parameters import (
    DomainParameters,
    RandomStrandGeneratorParameters,
    StrandDensityCrosslinkDistributerParameters,
)
from .networktype import NetworkType
from .network import Network
import numpy.random as npr
import numpy as np


def random_network(
    sizex,
    sizey,
    number_of_beads_per_strand,
    number_of_strands,
    contour_length_of_strand,
    crosslink_max_r,
    maximal_number_of_initial_crosslinks,
    crosslink_bin_size,
    seed=None,
    fix_boundary=False,
) -> Network:
    nt = NetworkType(
        DomainParameters(sizex, sizey, fix_boundary=fix_boundary),
        RandomStrandGenerator(
            RandomStrandGeneratorParameters(
                number_of_beads_per_strand=number_of_beads_per_strand,
                number_of_strands=number_of_strands,
                contour_length_of_strand=contour_length_of_strand,
            ),
            UniformStrandDistribution(sizex, sizey, seed),
        ),
        StrandDensityCrosslinkDistributer(
            StrandDensityCrosslinkDistributerParameters(
                crosslink_max_r,
                maximal_number_of_initial_crosslinks,
                number_of_beads_per_strand,
                number_of_strands,
                crosslink_bin_size,
            ),
            seed=seed,
        ),
        seed=seed,
    )
    return nt.generate()


def fibrin_network(
    sizex,
    sizey,
    number_of_beads_per_strand,
    number_of_strands,
    direction_spread,
    direction_angle,
    contour_length_of_strand,
    crosslink_max_r,
    maximal_number_of_initial_crosslinks,
    crosslink_bin_size,
    seed=None,
    fix_boundary=False,
    fix_boundary_north=False,
    fix_boundary_south=False,
    fix_boundary_east=False,
    fix_boundary_west=False,
) -> Network:
    nt = NetworkType(
        DomainParameters(
            sizex,
            sizey,
            fix_boundary=fix_boundary,
            fix_boundary_east=fix_boundary_east,
            fix_boundary_north=fix_boundary_north,
            fix_boundary_west=fix_boundary_west,
            fix_boundary_south=fix_boundary_south,
        ),
        RandomStrandGenerator(
            RandomStrandGeneratorParameters(
                number_of_beads_per_strand=number_of_beads_per_strand,
                number_of_strands=number_of_strands,
                contour_length_of_strand=contour_length_of_strand,
            ),
            VonMisesStrandDistribution(
                sizex, sizey, direction_angle, direction_spread, seed
            ),
        ),
        StrandDensityCrosslinkDistributerFast(
            StrandDensityCrosslinkDistributerParameters(
                crosslink_max_r,
                maximal_number_of_initial_crosslinks,
                number_of_beads_per_strand,
                number_of_strands,
                crosslink_bin_size,
            ),
            seed=seed,
        ),
        seed=seed,
        crosslink_angles=True,
    )
    return nt.generate()


def random_directed_network(
    sizex,
    sizey,
    number_of_beads_per_strand,
    number_of_strands,
    direction_spread,
    direction_angle,
    contour_length_of_strand,
    crosslink_max_r,
    maximal_number_of_initial_crosslinks,
    crosslink_bin_size,
    seed=None,
    fix_boundary=False,
    fix_boundary_north=False,
    fix_boundary_south=False,
    fix_boundary_east=False,
    fix_boundary_west=False,
) -> Network:
    nt = NetworkType(
        DomainParameters(
            sizex,
            sizey,
            fix_boundary=fix_boundary,
            fix_boundary_east=fix_boundary_east,
            fix_boundary_north=fix_boundary_north,
            fix_boundary_west=fix_boundary_west,
            fix_boundary_south=fix_boundary_south,
        ),
        RandomStrandGenerator(
            RandomStrandGeneratorParameters(
                number_of_beads_per_strand=number_of_beads_per_strand,
                number_of_strands=number_of_strands,
                contour_length_of_strand=contour_length_of_strand,
            ),
            VonMisesStrandDistribution(
                sizex, sizey, direction_angle, direction_spread, seed
            ),
        ),
        StrandDensityCrosslinkDistributer(
            StrandDensityCrosslinkDistributerParameters(
                crosslink_max_r,
                maximal_number_of_initial_crosslinks,
                number_of_beads_per_strand,
                number_of_strands,
                crosslink_bin_size,
            ),
            seed=seed,
        ),
        seed=seed,
    )
    return nt.generate()


def ISV_network(
    sizex,
    sizey,
    number_of_beads_per_strand,
    number_of_strands,
    contour_length_of_strand,
    crosslink_max_r,
    maximal_number_of_initial_crosslinks,
    crosslink_bin_size,
    spread_xaxis=1.0,
    seed=None,
    fix_boundary=None,
) -> Network:
    rng = np.random.default_rng(seed)
    if fix_boundary:
        domain = DomainParameters(sizex, sizey, fix_boundary=fix_boundary)
    else:
        domain = DomainParameters(
            sizex, sizey, fix_boundary_north=True, fix_boundary_south=True
        )
    nt = NetworkType(
        domain,
        RandomStrandGenerator(
            RandomStrandGeneratorParameters(
                number_of_beads_per_strand=number_of_beads_per_strand,
                number_of_strands=number_of_strands,
                contour_length_of_strand=contour_length_of_strand,
            ),
            StrandDistributionGeneral(
                lambda n: rng.normal(sizex * 0.5, spread_xaxis, n),
                lambda n: rng.uniform(0, sizey, n),
                lambda n: rng.uniform(0, 2 * np.pi, n),
            ),
        ),
        StrandDensityCrosslinkDistributer(
            StrandDensityCrosslinkDistributerParameters(
                crosslink_max_r,
                maximal_number_of_initial_crosslinks,
                number_of_beads_per_strand,
                number_of_strands,
                crosslink_bin_size,
            ),
            seed=seed,
        ),
        seed=seed,
    )
    return nt.generate()


def _random_laminin_positions(seed, beads_types, amount_of_laminin):
    rng = np.random.default_rng(seed)
    free_indices = [k for k, typ in enumerate(beads_types) if typ == "free"]
    list(rng.choice(free_indices, size=amount_of_laminin, replace=False))


from scipy.stats import truncnorm
from collections import defaultdict
import itertools

import matplotlib.pyplot as plt


def hexagonal(sizex, sizey, size):
    domain = DomainParameters(sizex, sizey, fix_boundary=True)

    horizontal_spacing = np.sqrt(3) * size  #  * (3.0 / 2.0)
    # vertical_spacing = size * (3.0/ 2.0)
    # vertical_spacing = size * (np.sqrt(3) / 2.0)
    vertical_spacing = size * (3.0 / 2.0)

    num_x = int(sizex / horizontal_spacing) + 1
    num_y = int(sizey / vertical_spacing) + 1

    coords = []
    ptypes = []
    c = np.sqrt(3) / 2
    s = 0.5

    bonds = []

    number_of_horizontal_beads = 2 * num_x + 1
    for r in range(0, num_y, 2):
        for i in range(2):
            for q in range(num_x):
                x = size * np.sqrt(3) * q + np.sqrt(3) * 0.5 * (r % 2)
                y = size * 3 * 0.5 * r
                index = len(coords)
                sign = -1 if i % 2 == 0 else 1

                if q == 0:
                    coords.extend(
                        [
                            (x - c * size, y + sign * s * size),  # 0
                        ]
                    )
                    ptypes.append("boundary")
                    bonds.append([index, index + 1])
                    index += 1
                else:
                    bonds.append([index - 1, index])

                coords.extend(
                    [
                        (x, y + sign * 1 * size),  # 1
                        (x + c * size, y + sign * s * size),  # 2
                    ]
                )
                if (
                    (r == 0 and i == 0)
                    or (r == num_y - 2 and i == 1)
                    or (r == num_y - 1 and i == 1)
                ):
                    ptypes.append("boundary")
                    ptypes.append("boundary")
                else:
                    ptypes.append("free")
                    if q == num_x - 1:
                        ptypes.append("boundary")
                    else:
                        ptypes.append("free")
                bonds.append([index, index + 1])
                # bonds.append([index+1, index+2])
            if i % 2 == 1:
                for q in range(0, num_x):
                    index = len(coords) - number_of_horizontal_beads + 2 * q
                    if q == 0:
                        bonds.append([index - number_of_horizontal_beads, index])
                    bonds.append([index - number_of_horizontal_beads + 2, index + 2])
            if r > 0 and i % 2 == 0:
                for q in range(0, num_x):
                    index = len(coords) - number_of_horizontal_beads + 2 * q + 1
                    bonds.append([index, index - number_of_horizontal_beads])
    assert len(ptypes) == len(coords)
    return Network(
        domain, coords, ptypes, bonds, ["polymer"] * len(bonds), [[0, 1, 2]], [0]
    )


def laminin(
    sizex,
    sizey,
    amount_of_laminin,
    network: Network,
    seed=None,
    x_dist_spread=1.0,
    y_dist_spread=0.0,
):

    pixel_to_bead = defaultdict(list)
    x_pos = []
    y_pos = []
    for k, pos in enumerate(network.beads_positions):
        nx, ny = int(pos[0]), int(pos[1])
        pixel_to_bead[nx, ny].append(k)
        x_pos.append(nx)
        y_pos.append(ny)

    locx = sizex // 2
    a_x = (0 - locx) / x_dist_spread
    b_x = (sizex - locx) / x_dist_spread
    x_pos = truncnorm.rvs(
        a_x, b_x, loc=locx, scale=x_dist_spread, size=amount_of_laminin
    )

    if y_dist_spread > 0:
        locy = sizey // 2
        a_y = (0 - locy) / y_dist_spread
        b_y = (sizey - locy) / y_dist_spread
        y_pos = truncnorm.rvs(
            a_y, b_y, loc=locy, scale=y_dist_spread, size=amount_of_laminin
        )
    else:
        y_pos = np.random.uniform(0, sizey, size=amount_of_laminin)

    # This can be smaller than amount_of_laminin
    free_beads_that_get_laminin_connection = itertools.chain(
        *[
            pixel_to_bead[int(x), int(y)]
            for (x, y) in zip(x_pos, y_pos)
            if (int(x), int(y)) in pixel_to_bead.keys()
        ]
    )

    # free_beads_that_get_laminin_connection = _random_laminin_positions(seed, network.beads_types, amount_of_laminin)
    laminin_positions = [
        network.beads_positions[k] for k in free_beads_that_get_laminin_connection
    ]
    laminin_types = ["boundary"] * len(laminin_positions)
    laminin_ids = [
        len(network.beads_positions) + k for k, _ in enumerate(laminin_positions)
    ]

    laminin_bonds = list(zip(free_beads_that_get_laminin_connection, laminin_ids))

    laminin_bonds_types = ["laminin"] * len(laminin_bonds)

    network.details_of_bondtypes["laminin"] = {
        "k": 1.0,  # ratio of spring_k,
        "r0": 0.0,  # used
    }

    network.beads_positions.extend(laminin_positions)
    network.beads_types.extend(laminin_types)

    network.bonds_groups.extend(laminin_bonds)
    network.bonds_types.extend(laminin_bonds_types)


def single_strand(
    sizex,
    sizey,
    start_x,
    start_y,
    angle,
    number_of_beads_per_strand,
    contour_length_of_strand,
    seed=None,
) -> Network:
    beads_to_middle_of_strand = number_of_beads_per_strand // 2
    length_single_bond = contour_length_of_strand / (number_of_beads_per_strand - 1)

    middle_of_strand_x = (
        start_x - np.cos(angle) * beads_to_middle_of_strand * length_single_bond
    )
    middle_of_strand_y = (
        start_y - np.sin(angle) * beads_to_middle_of_strand * length_single_bond
    )

    nt = NetworkType(
        DomainParameters(sizex, sizey, fix_boundary=True),
        RandomStrandGenerator(
            RandomStrandGeneratorParameters(
                number_of_beads_per_strand=number_of_beads_per_strand,
                number_of_strands=1,
                contour_length_of_strand=contour_length_of_strand,
            ),
            DeterministicStrandDistribution(
                [middle_of_strand_x], [middle_of_strand_y], [angle]
            ),
        ),
        None,
        seed=seed,
    )
    return nt.generate()


from .regular_network import (
    RegularNetwork,
    RegularNetworkParameters,
    RegularCrosslinker,
)


def regular(
    sizex,
    sizey,
    number_of_fibers_per_side,
    number_of_beads_per_strand,
    fix_boundary,
    single_side=False,
):
    print("single_side = ", single_side)
    par = RegularNetworkParameters(
        number_of_fibers_per_side,
        number_of_beads_per_strand,
        only_vertical_strands=single_side,
    )

    cross = RegularCrosslinker(par) if not single_side else None

    nt = NetworkType(
        DomainParameters(sizex, sizey, fix_boundary),
        RegularNetwork(par),
        cross,
    )

    return nt.generate()


def single_spring(
    sizex,
    sizey,
    number_of_strands,
    number_of_beads_per_strand,
    contour_length_of_strand,
    seed=None,
) -> Network:

    assert number_of_strands % 2 == 1
    n = (number_of_strands - 1) // 2
    x_pos = [100 + 100.0 * k / (2.0 * n) for k in range(number_of_strands)]
    angle_up = np.arctan(np.sqrt(contour_length_of_strand**2 - (1 / n) ** 2) * 2 * n)
    angle_down = np.arctan(-np.sqrt(contour_length_of_strand**2 - (1 / n) ** 2) * 2 * n)
    angles = [angle_up, angle_down] * ((number_of_strands - 1) // 2) + [angle_up]

    nt = NetworkType(
        DomainParameters(sizex, sizey, fix_boundary=True),
        RandomStrandGenerator(
            RandomStrandGeneratorParameters(
                number_of_beads_per_strand=number_of_beads_per_strand,
                number_of_strands=number_of_strands,
                contour_length_of_strand=contour_length_of_strand,
            ),
            DeterministicStrandDistribution(x_pos, [100] * number_of_strands, angles),
        ),
        TipToTailCrosslinkDistributer(
            number_of_beads_per_strand=number_of_beads_per_strand,
            number_of_strands=number_of_strands,
        ),
        seed=seed,
    )
    return nt.generate()


def two_crosslinked_strands(r0):
    beads = 4
    strand1 = single_strand(10, 10, 5, 5, 0, beads, (beads - 1) * r0)
    strand2 = single_strand(10, 10, 5, 4, 0, beads, (beads - 1) * r0)
    net = strand1 + strand2

    crosslinker = DeterministicCrosslinkDistributer([[1, 5]])
    crosslinker.distribute_crosslinkers(net)
    crosslinker.add_crosslink_angles(net)

    print(net)

    return net
