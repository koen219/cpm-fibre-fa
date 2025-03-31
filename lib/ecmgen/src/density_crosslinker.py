from .crosslink_distributors import CrosslinkDistributer
from .network import BOND, BONDTYPE, BEADID, FIBREID, BONDID, Network
from .binner import FiberBin
from .parameters import StrandDensityCrosslinkDistributerParameters
from scipy.signal import convolve2d
from scipy.special import binom

import numpy as np
from typing import Optional, List, Tuple, Dict, Set
from abc import ABC, abstractmethod
import math
import itertools
from collections import defaultdict
from .crosslink_distributors import _CrosslinkQuantizer


class StrandDensityCrosslinkDistributer(CrosslinkDistributer):
    def __init__(
        self, par: StrandDensityCrosslinkDistributerParameters, seed: Optional[int]
    ):
        self._par = par
        self._rng = np.random.default_rng(seed)
        self._binner: FiberBin

        self._quantizer = _CrosslinkQuantizer(self._par.crosslink_max_r, 10)

        # self.spring_options = self._quantizer.spring_options(self.par.crosslink_k)
        # {'crosslinker': dict(
        #    k=par.crosslink_k, r0=par.crosslink_r0)}

    def select_bonds(self, network: Network) -> List[Tuple[BOND, BONDTYPE]]:
        pos = network.beads_positions

        density_bin = self._makeDensityBin(network)
        pairings = self._pairings(density_bin)
        total_num_pairings = np.sum(pairings)

        # for each sampled bin, generate the appropriate number of crosslinks by drawing without replacement, with uniform probability
        sample = self._sample_from_2d_array(
            pairings / total_num_pairings,
            size=self._par.maximal_number_of_initial_crosslinks,
        )

        selected_bonds = self._select_bonds_per_fibre(network, sample)
        selected_types = []

        for bond in selected_bonds:
            r = np.linalg.norm(np.array(pos[bond[0]]) - np.array(pos[bond[1]]))
            typ = self._quantizer.computetype(float(r))
            selected_types.append(typ)

            network.details_of_bondtypes[typ] = {
                "r0": self._quantizer.spring_options(0.0)[typ]["r0"],
                "k": 1,
            }

        return list(zip(selected_bonds, selected_types))

    def _map_bead_to_fiber(self, bead: BEADID) -> FIBREID:
        b = self._par.number_of_beads_per_strand
        return math.floor(bead / b)

    def _makeDensityBin(self, network: Network):
        # Function calculating the local densities which are used for the calculation of crosslinking probabilty
        # bin the network - note that the binning lattice is not the same as the CPM lattice!
        bin_size_x = self._par.crosslink_bin_size
        bin_size_y = self._par.crosslink_bin_size
        num_bins_x = int(network.domain.sizex / bin_size_x)
        num_bins_y = int(network.domain.sizey / bin_size_y)

        crosslink_binner = FiberBin(
            num_bins_x,
            num_bins_y,
            bin_size_x,
            bin_size_y,
            self._par.number_of_beads_per_strand,
            self._par.number_of_strands,
            network.domain.sizex,
            network.domain.sizey,
        )

        bondgroup = np.array(network.bonds_groups)
        particlepos = np.array(network.beads_positions)
        crosslink_binner.bin_network(bondgroup, particlepos)
        density_bin = np.array(crosslink_binner.density_bin)

        self._binner = crosslink_binner
        return density_bin

    def _pairings(self, density_bin):
        """
        Based on local densities, count the number of possible pairings (same fiber pairing is included)
        """
        filt = np.ones((3, 3))  # Moore NBH
        summ = convolve2d(density_bin, filt)
        # this is n choose 2 with n = total number of bonds in the moore nbh
        pairings = 0.5 * (summ - 1) * summ
        return pairings

    def _find_bonds_in_neighbourhood(self, nx, ny) -> Set[BONDID]:
        """
        Returns the bondids of all bonds in a Moore neighbourhood of nx,ny
        """
        bonds_in_nbhd: Set[BONDID] = set()
        for dx, dy in itertools.product([-1, 0, 1], [-1, 0, 1]):
            try:
                bondstoadd = self._binner.parent_bin[ny + dy][nx + dx]
            except IndexError as err:  # This takes care of the edges of the bin
                bondstoadd = {}
            bonds_in_nbhd = bonds_in_nbhd.union(bondstoadd)
        return bonds_in_nbhd

    def _find_all_bead_pairs_with_different_fibers(
        self, network: Network, bonds_in_nbhd: Set[BONDID]
    ):
        """
        Pairs all beads in the list of bonds that are connected to different fibers.
        """
        # We create a list of vertex_pairs such that no two vertices are part of the same fibre.
        fiber_dict: Dict[FIBREID, List[BEADID]] = defaultdict(list)

        for bondid in bonds_in_nbhd:
            for bead in network.bonds_groups[bondid]:
                fiber_id = self._map_bead_to_fiber(bead)
                fiber_dict[fiber_id].append(bead)

        # Example: fiber_dict = { 1: [1,2,3], 2: [4,5,6]: 3: [7,8] } then combinations(fiber_dict.values(),2) =
        # [ ([1,2,3] ,[4,5,6]), ([1,2,3], [7,8]), ([4,5,6], [7,8])]
        # So fiberpair = tuple of list of beads that are guaranteed to be on different fibers

        all_bead_pairs = []
        fiberpair: Tuple[List[BEADID], List[BEADID]]
        for fiberpair in itertools.combinations(fiber_dict.values(), 2):
            beadpair: List[BOND] = list(itertools.product(*fiberpair))
            all_bead_pairs.extend(beadpair)
        return all_bead_pairs

    def _crosslink_r_dist(self, r):  # distribution of possible radii of a crosslinker
        if r <= self._par.crosslink_max_r:
            p = 2 / self._par.crosslink_max_r - 2 / self._par.crosslink_max_r**2 * r
        else:
            p = 0
        return p

    def _sample_bonds_on_distance(self, bead_pairs, num_of_samples, network: Network):
        dist = []
        pos = np.array(network.beads_positions)
        cum_p = 0
        for k0, k1 in bead_pairs:
            r = np.linalg.norm(np.array(pos[k0]) - np.array(pos[k1]))  # type: ignore
            p = self._crosslink_r_dist(r)
            dist.append(p)
            cum_p += p
        if cum_p == 0:  # Could not make a crosslinker
            return []
        # normalize
        prob = [p / cum_p for p in dist]

        # sample the vertex pairs -- closer pairs are more likely to be chosen
        pairs_indices = self._rng.choice(
            len(bead_pairs),
            size=num_of_samples,
            replace=False,
            p=prob,
        )
        return [
            sorted(bead_pairs[k]) for k in pairs_indices  # We should have used sets
        ]

    def _select_bonds_per_fibre(self, network: Network, sample):
        # For each selected site, try to create a crosslinker that is not connected to the same fiber.
        selected_bonds = {}
        allready_crosslinked_beads = list()
        for ny, nx in sample:
            # Consider the bins next to the current bin and take all bonds that are in there as well.
            bonds_in_nbhd = self._find_bonds_in_neighbourhood(nx, ny)
            all_bead_pairs = self._find_all_bead_pairs_with_different_fibers(
                network, bonds_in_nbhd
            )
            # remove beads that already have a crosslinker
            filtered_bead_pairs = list(
                filter(
                    lambda bp: not (
                        bp[0] in allready_crosslinked_beads
                        or bp[1] in allready_crosslinked_beads
                    ),
                    all_bead_pairs,
                )
            )

            if len(filtered_bead_pairs) == 0:  # If there was only one fiber
                continue
            if (
                len(
                    bonds := self._sample_bonds_on_distance(
                        filtered_bead_pairs, 1, network
                    )
                )
                > 0
            ):
                allready_crosslinked_beads.extend(itertools.chain(*bonds))
                selected_bonds[(ny, nx)] = bonds[0]
        return list(selected_bonds.values())

    def _sample_from_2d_array(self, array, size):
        """
        Sample 'size' indices from array with weights as in the array.
        """
        # Create a flat copy of the array
        flat = array.flatten()

        # Then, sample an index from the 1D array with the
        # probability distribution from the original array
        sample_index = np.random.choice(a=flat.size, p=flat, replace=True, size=size)

        # Take this index and adjust it so it matches the original array
        adjusted_index = np.unravel_index(sample_index, array.shape)
        return list(zip(*adjusted_index))  # type: ignore


class StrandDensityCrosslinkDistributerFast(CrosslinkDistributer):
    def __init__(
        self, par: StrandDensityCrosslinkDistributerParameters, seed: Optional[int]
    ):
        self._par = par
        self._rng = np.random.default_rng(seed)
        self._binner: FiberBin

        self._quantizer = _CrosslinkQuantizer(self._par.crosslink_max_r, 10)

        # self.spring_options = self._quantizer.spring_options(self.par.crosslink_k)
        # {'crosslinker': dict(
        #    k=par.crosslink_k, r0=par.crosslink_r0)}

    def same_fiber(self, bead1, bead2) -> bool:
        b = self._par.number_of_beads_per_strand
        return math.floor(bead1 / b) == math.floor(bead2 / b)

    def select_bonds(self, network: Network) -> List[Tuple[BOND, BONDTYPE]]:
        bin = self._binNetwork(network)
        pos = np.asarray(network.beads_positions)

        bonds, types = list(), list()

        number_of_combinations = sum(
            binom(len(ids), 2) for ids in bin.values() if len(ids) >= 2
        )
        if number_of_combinations == 0:
            return []
        prob = self._par.maximal_number_of_initial_crosslinks / number_of_combinations
        print(
            self._par.maximal_number_of_initial_crosslinks, number_of_combinations, prob
        )
        for ids in bin.values():
            if len(ids) < 2:
                continue
            for p1, p2 in itertools.combinations(ids, 2):
                if self.same_fiber(p1, p2):
                    continue
                if self._rng.random() > prob:
                    continue
                new_bond = (p1, p2)
                r = np.linalg.norm(pos[p1] - pos[p2])
                bond_type = self._quantizer.computetype(float(r))

                bonds.append(new_bond)
                types.append(bond_type)

        for typ in set(types):
            network.details_of_bondtypes[typ] = {
                "r0": self._quantizer.spring_options(0.0)[typ]["r0"],
                "k": 1,
            }

        return list(zip(bonds, types))

    def _binNetwork(self, network: Network):
        """
        Bins 2D positions into a dictionary.

        Args:
            positions (numpy.ndarray): An Nx2 array of positions.
            bin_size (float): Size of each bin.
            domain_size (tuple): Maximum x, y, z coordinates defining the domain.

        Returns:
            dict: A dictionary where keys are bin indices (tuples) and values are lists of position indices.
        """
        positions = np.asarray(network.beads_positions)
        domain_size = (network.domain.sizex, network.domain.sizey)
        bin_size = self._par.crosslink_bin_size

        # Filter out positions outside the domain_size
        within_domain = np.all((positions >= 0) & (positions < domain_size), axis=1)
        positions = positions[within_domain]
        valid_indices = np.where(within_domain)[0]

        # Compute bin indices for each valid position
        bin_indices = (positions // bin_size).astype(int)

        # Group positions by bins using a dictionary
        bins = defaultdict(list)
        for idx, bin_idx in zip(valid_indices, bin_indices):
            bins[tuple(bin_idx)].append(idx)

        return bins
