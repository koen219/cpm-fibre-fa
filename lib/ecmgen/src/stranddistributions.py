import numpy.typing as npt
import numpy as np
from abc import ABC, abstractmethod
from typing import Optional, Callable


class StrandDistribution(ABC):
    @abstractmethod
    def pos_x_dist(self, n) -> npt.NDArray[np.float64]:
        pass

    @abstractmethod
    def pos_y_dist(self, n) -> npt.NDArray[np.float64]:
        pass

    @abstractmethod
    def angle_dist(self, n) -> npt.NDArray[np.float64]:
        pass

class StrandDistributionGeneral(StrandDistribution):
    def __init__(
            self,
    pos_x_dist: Callable[[int], npt.NDArray[np.float64]],
    pos_y_dist: Callable[[int], npt.NDArray[np.float64]],
    angle_dist: Callable[[int], npt.NDArray[np.float64]],
    ):
        self._pos_x_dist = pos_x_dist
        self._pos_y_dist = pos_y_dist
        self._angle_dist = angle_dist

    def pos_x_dist(self, n) -> npt.NDArray[np.float64]:
        return self._pos_x_dist(n)

    def pos_y_dist(self, n) -> npt.NDArray[np.float64]:
        return self._pos_y_dist(n)

    def angle_dist(self, n) -> npt.NDArray[np.float64]:
        return self._angle_dist(n)


class UniformStrandDistribution(StrandDistribution):
    def __init__(self, sizex, sizey, seed: Optional[int] = None):
        self._sizex = sizex
        self._sizey = sizey
        if seed:
            self._rng = np.random.default_rng(seed)
        else:
            self._rng = np.random.default_rng()

    def pos_x_dist(self, n):
        """
        Return starting positions of beads
        """
        return self._rng.uniform(low=0, high=self._sizex, size=n)

    def pos_y_dist(self, n):
        return self._rng.uniform(0, self._sizey, size=n)

    def angle_dist(self, n):
        return self._rng.uniform(0, 2 * np.pi, size=n)

class VonMisesStrandDistribution(StrandDistribution):
    def __init__(self, sizex, sizey, mu, kappa, seed= None):
        self._uniform_strand_distribution = UniformStrandDistribution(sizex,sizey,seed)
        self._mu = mu
        self._kappa = kappa
        if seed:
            self._rng = np.random.default_rng(seed)
        else:
            self._rng = np.random.default_rng()
    
    def pos_x_dist(self, n):
        return self._uniform_strand_distribution.pos_x_dist(n)

    def pos_y_dist(self, n):
        return self._uniform_strand_distribution.pos_y_dist(n)

    def angle_dist(self, n):
        return self._rng.vonmises(self._mu, self._kappa, size=n)

class DeterministicStrandDistribution(StrandDistribution):
    def __init__(self, x_samples, y_samples, angle_samples):
        self._x_samples = np.array(x_samples)
        self._y_samples = np.array(y_samples)
        self._angle_samples = np.array(angle_samples)

    def pos_x_dist(self, n):
        """
        Return starting positions of beads
        """
        if len(self._x_samples) != n:
            raise ValueError(
                f"There are {n} samples requested but we only have {len(self._x_samples)}"
            )
        return self._x_samples

    def pos_y_dist(self, n):
        if len(self._y_samples) != n:
            raise ValueError(
                f"There are {n} samples requested but we only have {len(self._y_samples)}"
            )
        return self._y_samples

    def angle_dist(self, n):
        if len(self._angle_samples) != n:
            raise ValueError(
                f"There are {n} samples requested but we onlangle have {len(self._angle_samples)}"
            )
        return self._angle_samples
