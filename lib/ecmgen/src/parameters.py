from dataclasses import dataclass, field


@dataclass
class DomainParameters:
    sizex: int
    sizey: int

    fix_boundary: bool = field(default=False)
    fix_boundary_north: bool = field(default=False)
    fix_boundary_south: bool = field(default=False)
    fix_boundary_west: bool = field(default=False)
    fix_boundary_east: bool = field(default=False)
    @property
    def Lx(self):
        return self.sizex // 2

    @property
    def Ly(self):
        return self.sizey // 2


@dataclass
class RandomStrandGeneratorParameters:
    number_of_beads_per_strand: int
    number_of_strands: int
    contour_length_of_strand: float


@dataclass
class StrandDensityCrosslinkDistributerParameters:
    crosslink_max_r: float
    maximal_number_of_initial_crosslinks: int
    number_of_beads_per_strand: int
    number_of_strands: int

    crosslink_bin_size: float
