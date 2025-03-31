from abc import ABC, abstractmethod
from typing import Optional, List, Tuple, Dict, Set
from .network import Network
from .network import BOND, BONDTYPE, BEADID, FIBREID, BONDID, Network
import itertools
import numpy as np

import logging

_logger = logging.getLogger(__name__)


class _CrosslinkQuantizer:
    def __init__(self, crosslink_max_r, number):
        self._max = crosslink_max_r
        self._num = number

        quantizations = np.linspace(0, self._max, self._num)

        self.types = [f"cross_{i}" for i in range(self._num)]
        self.types_r0 = [q for q in quantizations]

    def computetype(self, R: float):
        if R > self._max:
            return self.types[-1]
        if R < 0:
            return 0

        r = round((R / self._max) * self._num)

        return self.types[r - 1]

    def spring_options(self, stiffness):
        return {
            typ: dict(r0=r0, k=stiffness) for typ, r0 in zip(self.types, self.types_r0)
        }


class CrosslinkDistributer(ABC):
    def distribute_crosslinkers(self, network: Network):
        bonds_to_add = list()
        types_to_add = list()

        beads_with_crosslinker = set()

        selected_bonds_and_types = self._bonds_to_crosslink(network)

        for bond, bond_typ in selected_bonds_and_types:
            sorted_bond: List[int] = sorted(bond)

            if (
                sorted_bond in bonds_to_add
                or sorted_bond[0] in beads_with_crosslinker
                or sorted_bond[1] in beads_with_crosslinker
            ):
                continue

            bonds_to_add.append(bond)
            types_to_add.append(bond_typ)
            beads_with_crosslinker.add(bond[0])
            beads_with_crosslinker.add(bond[1])

        network.bonds_groups.extend(bonds_to_add)
        network.bonds_types.extend(types_to_add)

    def add_crosslink_angles(self, network: Network):
        print("Add crosslink angles")
        selected_bonds_and_types = self._bonds_to_crosslink(network)

        # For each bond, we have to add an angle constraint to some triples of neighbours.
        # I only add an additional angle if the crosslink is in the middle of a fiber.
        # The ends of fibers remain free to rotate. Partly because it is unclear to which
        # end a angle should be added
        angles_to_add = []
        thetas_to_add = []

        beads_positions = np.array(network.beads_positions)
        bonds_group = dict()
        for bond, bondtype in zip(network.bonds_groups, network.bonds_types):
            if bondtype == "polymer":
                bonds_group.setdefault(bond[0], []).append(bond[1])
                bonds_group.setdefault(bond[1], []).append(bond[0])
        print(f"{bonds_group=}")

        for bond, bondtype in selected_bonds_and_types:
            angles = [
                (bond[0], bond[1], bonds_group.get(bond[1])[0]),
                (bond[1], bond[0], bonds_group.get(bond[0])[0]),
            ]

            print(list(angles))
            for a, b, c in angles:
                # a, b, c = next(angles)
                # a, b, c = np.random.choice(angles, 1)
                # for a, b, c in angles:
                x_b = beads_positions[b]
                v_a = beads_positions[a] - x_b
                v_c = beads_positions[c] - x_b

                norm_v_a = np.linalg.norm(v_a)
                norm_v_c = np.linalg.norm(v_c)

                if norm_v_a == 0 or norm_v_c == 0:
                    continue

                if norm_v_a < 1.0 or norm_v_c == 1.0:
                    continue
                argument = np.dot(v_a, v_c) / (norm_v_a * norm_v_c)
                if argument < -1.0:
                    argument = -1.0
                if argument > 1.0:
                    argument = 1.0

                theta = np.arccos(np.dot(v_a, v_c) / (norm_v_a * norm_v_c))
                angle = (a, b, c)

                if theta > np.pi * 0.5:
                    continue

                angles_to_add.append(angle)
                thetas_to_add.append(theta)

        quantizer = _CrosslinkQuantizer(np.pi * 0.5, 20)
        angle_types = list()
        for theta in thetas_to_add:
            typ_name = quantizer.computetype(theta)
            t0 = quantizer.spring_options(0)[typ_name]["r0"]
            print(f"{theta=} -> type={typ_name} with {t0=}")
            angle_types.append("angle_" + typ_name)
            network.details_of_angletypes["angle_" + typ_name] = {
                "t0": t0,
                "k": 1,
            }
        print(network.details_of_angletypes)
        print(f"Adding {len(angles_to_add)}")
        network.angle_groups.extend(angles_to_add)
        network.angle_types.extend(angle_types)

    def _bonds_to_crosslink(self, network: Network):
        if not hasattr(self, "_selected_bonds_and_types"):
            self._selected_bonds_and_types = self.select_bonds(network)
        return self._selected_bonds_and_types

    @abstractmethod
    def select_bonds(self, network: Network) -> List[Tuple[BOND, BONDTYPE]]:
        pass


class DeterministicCrosslinkDistributer(CrosslinkDistributer):
    def __init__(self, crosslinks: List[BOND]):
        self._crosslinks = crosslinks

    def select_bonds(self, network: Network):
        types = ["cross"] * len(self._crosslinks)
        network.details_of_bondtypes["cross"] = {"r0": 1, "k": 1}

        return list(zip(self._crosslinks, types))


class TipToTailCrosslinkDistributer(CrosslinkDistributer):
    def __init__(self, number_of_beads_per_strand, number_of_strands):
        self._num_beads_per_strand = number_of_beads_per_strand
        self._num_strands = number_of_strands

    def select_bonds(self, network: Network):
        selected_bonds = []
        crosslink_type = "crosslinker"

        network.details_of_bondtypes[crosslink_type] = {"r0": 0, "k": 1}

        for i in range(self._num_strands - 1):
            b0 = self._num_beads_per_strand * i + self._num_beads_per_strand - 1
            b1 = b0 + 1
            selected_bonds.append(((b0, b1), crosslink_type))

        return selected_bonds
