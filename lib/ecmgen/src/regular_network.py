from .strandgens import StrandGenerator
from .crosslink_distributors import CrosslinkDistributer
from .network import Network, DomainParameters

import numpy as np
from dataclasses import dataclass,field

@dataclass
class RegularNetworkParameters:
    number_of_strands: int
    number_of_beads_per_strand: int
    
    only_vertical_strands: bool = field(default=False)

class RegularNetwork(StrandGenerator):
    def __init__(self, par: RegularNetworkParameters):
        self._par = par
        self._network = None

    def fix_boundaries(self, network: Network):
        typeid = network.beads_types

        for k,_ in enumerate(typeid):
            bead = k % self._par.number_of_beads_per_strand
            if bead == 0 or bead == self._par.number_of_beads_per_strand -1:
                network.beads_types[k] = "boundary"

    def build_strands(self, network: Network) -> Network:
        ###############################
        #  Generate network: strands  #
        ###############################

        # init_network_data is data that is needed to initalize the hoomd snapshot
        # it is deleted in set_hoomd_wrap. After that we save all changes in the data directly to the snapshot
        #

        pos, typeid = self._pos_gen(network.domain)
        # pos = self._rotate_network(pos)
        bondsgroup, bondstype = self._bond_gen(network.domain)
        anglegroup, angletypes = self._angle_gen(network.domain)

        network.beads_positions.extend(pos)
        network.beads_types.extend(typeid)

        network.bonds_groups.extend(bondsgroup)
        network.bonds_types.extend(bondstype)

        network.angle_groups.extend(anglegroup)
        network.angle_types.extend(angletypes)

        return network


    def _vertical_coordinate(self, strand, bead, box_x, box_y):
        Lx = 0 #box_x / 2
        Ly = 0 #box_y / 2
        strands = self._par.number_of_strands
        beads = self._par.number_of_beads_per_strand
        strands_per_direction = strands if self._par.only_vertical_strands else strands // 2
        x = (strand / strands_per_direction) * box_x - Lx
        y = (bead / beads) * box_y - Ly
        return x, y

    def _horizonal_coordinate(self, strand, bead, box_x, box_y):
        Lx = 0 # box_x / 2
        Ly = 0 # box_y / 2
        strands = self._par.number_of_strands
        beads = self._par.number_of_beads_per_strand
        strands_per_direction = strands // 2
        x = (bead / beads) * box_x - Lx
        y = ((strand - strands_per_direction) /
             strands_per_direction)*box_y - Ly
        return x, y

    def _pos_gen(self, domain: DomainParameters):
        num_particles = self._par.number_of_strands*self._par.number_of_beads_per_strand
        pos = [] # np.zeros((num_particles, 3))
        typeid = []
        
        if self._par.only_vertical_strands:
            strands_per_direction = self._par.number_of_strands
        else:
            strands_per_direction = self._par.number_of_strands//2

        for strand in range(self._par.number_of_strands):
            for bead in range(self._par.number_of_beads_per_strand):
                index = strand*self._par.number_of_beads_per_strand + bead

                if strand < strands_per_direction:  # first half vertical
                    x, y = self._vertical_coordinate(strand, bead, domain.sizex, domain.sizey)
                else: # if only vertical strands, this part does not get activated
                    x, y = self._horizonal_coordinate(strand, bead, domain.sizex, domain.sizey)

                if bead == 0 or bead == self._par.number_of_beads_per_strand - 1:  # begin and endpoints are fixed
                    typ = 'boundary'
                else:
                    typ = 'free'
                pos.append([x,y])
                typeid.append(typ)

        return pos, typeid

    def _bond_gen(self, domain):
        beads = self._par.number_of_beads_per_strand
        strands = self._par.number_of_strands
        bondsgroup = np.empty(
            shape=(strands*(beads-1), 2), dtype=int)
        bonds_indices = np.repeat(
            beads*np.arange(0, strands, dtype=int), beads-1)
        bondsgroup[:, 0] = bonds_indices + \
            np.tile(np.arange(0, beads-1,
                              1, dtype=int), strands)
        bondsgroup[:, 1] = bonds_indices + \
            np.tile(np.arange(1, beads, 1,
                    dtype=int), strands)
        bondsgroup = bondsgroup.tolist()
        bondstype = ['polymer'] * len(bondsgroup)
        return bondsgroup, bondstype

    # Angle generation
    def _angle_gen(self, domain):
        beads = self._par.number_of_beads_per_strand
        strands = self._par.number_of_strands

        angles_group = np.empty(
            shape=(strands*(beads-2), 3), dtype=int)
        angles_indices = np.repeat(
            beads*np.arange(0, strands, dtype=int), beads-2)
        angles_group[:, 0] = angles_indices + \
            np.tile(np.arange(0, beads-2,
                              1, dtype=int), strands)
        angles_group[:, 1] = angles_indices + \
            np.tile(np.arange(1, beads-1,
                              1, dtype=int), strands)
        angles_group[:, 2] = angles_indices + \
            np.tile(np.arange(2, beads, 1,
                    dtype=int), strands)
        anglegroup = angles_group.tolist()
        angles_types = [0] * len(anglegroup)
        return anglegroup, angles_types

class RegularCrosslinker(CrosslinkDistributer):
    """
    Creates crosslinkers on a square grid.
    """

    def __init__(self, par: RegularNetworkParameters):
        self._num_beads = par.number_of_beads_per_strand
        self._num_strands = par.number_of_strands

    def select_bonds(self, network: Network):
        selected_bonds = []
        crosslink_type = "crosslinker"

        network.details_of_bondtypes[crosslink_type] = {'r0': 0, 'k': 1}

        assert self._num_beads == (
            self._num_strands // 2), "Can not crosslink non square grid"
        d = self._num_beads

        selected_bonds = list(zip(
            np.arange(d*d, dtype=int).reshape((d, d)).transpose().flatten(),
            d*d + np.arange(d*d, dtype=int)
        ))

        return [ (bond, crosslink_type) for bond in selected_bonds ]