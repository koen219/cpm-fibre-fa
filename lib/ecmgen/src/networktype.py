import numpy as np
import numpy.typing as npt
from typing import Optional
from .strandgens import StrandGenerator
from .parameters import DomainParameters
from .crosslink_distributors import CrosslinkDistributer

from .network import Network

import logging

logger = logging.getLogger(__name__)


class NetworkType:
    """Config class of the network. Used to generate a Network class."""

    def __init__(
        self,
        domain: DomainParameters,
        strandgenerator: StrandGenerator,
        crosslink_distributor: Optional[CrosslinkDistributer],
        seed: Optional[int] = None,
        crosslink_angles=False,
    ):
        self._strand_generator: StrandGenerator = strandgenerator
        self._crosslink_distributor = crosslink_distributor
        self._rng = np.random.default_rng(seed=seed)
        self._network = Network(domain)
        self._crosslink_angles = crosslink_angles

        logger.info(
            "Initiate NetworkType with %s and crosslinked with %s"
            % (strandgenerator, crosslink_distributor)
        )

    def generate(self) -> Network:
        """Generates a network from the generators. Throws exceptions when some network are not neatly generated."""
        self._strand_generator.build_strands(self._network)
        # if self._network.domain.fix_boundary:
        self._strand_generator.fix_boundaries(self._network)
        if self._crosslink_distributor:
            self._crosslink_distributor.distribute_crosslinkers(self._network)
            if self._crosslink_angles:
                self._crosslink_distributor.add_crosslink_angles(self._network)
        return self._network

    @property
    def network(self):
        return self._network
