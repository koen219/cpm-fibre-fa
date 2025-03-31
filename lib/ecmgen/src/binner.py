import numpy as np
import math
import operator
import copy
from scipy.stats import binned_statistic_2d
from itertools import product


class FiberBin:
    def __init__(
        self,
        num_bins_x,
        num_bins_y,
        bin_size_x,
        bin_size_y,
        beads,
        strands,
        sizex,
        sizey,
        offset_x=0,
        offset_y=0,
        precision=8,
    ):
        self.num_bins_x = num_bins_x
        self.num_bins_y = num_bins_y
        self.bin_size_x = bin_size_x
        self.bin_size_y = bin_size_y
        self.sizex = sizex
        self.sizey = sizey
        self.offset_x = offset_x
        self.offset_y = offset_y
        self.precision = precision
        self.num_beads = beads
        self.num_strands = strands
        # initialize default coordinate conversion
        self.set_coord_to_bin()

        ## initialize the bins
        # self.initialize_bins()

    def initialize_bins(self):
        # self.segment_bin = self.generate_bin([]) Not used
        self.parent_bin = self.generate_bin([])
        self.density_bin = self.generate_bin(0)  # Not used

    def generate_bin(self, inner_obj=[]):
        # build empy bin [[inner_obj]] and return it
        newbin = []
        for ny in range(self.num_bins_y):
            tmplist = []
            for nx in range(self.num_bins_x):
                tmplist.append(copy.copy(inner_obj))
            newbin.append(tmplist)

        return newbin

    def set_coord_to_bin(self, f=None):
        # define custom conversion from Euclidean coordinates to bin index
        if f is None:
            self._coord_to_bin_f = self._default_coord_to_bin
        else:
            self._coord_to_bin_f = f

    def _default_coord_to_bin(self, v):
        # default coord_to_fin method
        # takes coordinates w=(x,y) and returns bin index (nx,ny)
        x, y = v
        nx = int(x / self.bin_size_x)
        ny = int(y / self.bin_size_y)
        return nx, ny

    def coord_to_bin(self, v):
        # takes coordinates v=(x,y) and returns bin index (nx,ny)
        # raises exception if coordinates does not fall within one of the bins
        nx, ny = self._coord_to_bin_f(v)
        if nx < 0 or ny < 0 or nx >= self.num_bins_x or ny >= self.num_bins_y:
            raise ValueError
        return nx, ny

    def _interpolate_beads(self, pos, interpolation_number):
        inter = np.zeros(
            (interpolation_number * self.num_strands * (self.num_beads - 1), 2)
        )
        t = np.linspace(0, 1, interpolation_number, endpoint=False)
        for strand in range(self.num_strands):
            for k in range(self.num_beads)[:-1]:
                bondid = strand * (self.num_beads - 1) + k
                beadid = strand * self.num_beads + k
                P = pos[beadid]
                Q = pos[beadid + 1]
                inter[
                    bondid * interpolation_number : (bondid + 1) * interpolation_number,
                    0,
                ] = P[0] + t * (Q[0] - P[0])
                inter[
                    bondid * interpolation_number : (bondid + 1) * interpolation_number,
                    1,
                ] = P[1] + t * (Q[1] - P[1])
        return inter

    def bin_network(self, bonds_group, pos):
        # bins the network

        # initialize / clear the bins
        self.initialize_bins()
        N = 100
        inter = self._interpolate_beads(pos, N)
        X_inter = inter[:, 0]
        Y_inter = inter[:, 1]
        bin_x = np.linspace(
            0.0 - self.offset_x, self.sizex - self.offset_x, self.num_bins_x
        )
        bin_y = np.linspace(
            0.0 - self.offset_y, self.sizey - self.offset_y, self.num_bins_y
        )
        num_bin_x = len(bin_x) - 1
        num_bin_y = len(bin_y) - 1
        binning_range = np.array(
            [
                [0.0 - self.offset_x, self.sizex - self.offset_x],
                [0.0 - self.offset_y, self.sizey - self.offset_y],
            ]
        )
        counts, x_edge, y_edge, binnumber = binned_statistic_2d(
            X_inter,
            Y_inter,
            values=[],
            statistic="count",
            bins=[bin_x, bin_y],
            expand_binnumbers=True,
            range=binning_range,
        )
        self.counts = counts
        for k, [nx, ny] in enumerate(binnumber.T):
            if nx >= self.num_bins_x or ny >= self.num_bins_y:
                self.parent_bin[ny - 1][nx - 1] = []
                continue
            #            if nx == self.num_bins_x or ny == self.num_bins_y:
            #                continue
            #            self.density_bin[ny-1][nx-1] = counts[nx-1][ny-1]
            #            self.parent_bin[ny-1][nx-1].append(k// N)
            parent_bin_set = set(self.parent_bin[ny - 1][nx - 1])
            parent_bin_set.add(k // N)
            self.parent_bin[ny - 1][nx - 1] = list(parent_bin_set)
        for nx, ny in product(range(self.num_bins_x), range(self.num_bins_y)):
            self.density_bin[ny][nx] = len(self.parent_bin[ny][nx])
