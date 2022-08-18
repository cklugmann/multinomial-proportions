import warnings
from typing import Union, Optional, List, Callable, Any

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


def requires_fit(fn: Callable) -> Callable:
    def closure(self, *args, **kwargs) -> Any:
        if "sample_size" not in self.observed:
            raise ValueError("You need to call the `fit()` method first!")

        return fn(self, *args, **kwargs)

    return closure


class MultinomialEstimator:

    """
    This class can be used to calculate simultaneous confidence intervals for multinomial proportions or to determine
    the sample size in order to specify the confidence intervals for a given precision.
    """

    def __init__(self, num_cells: int, tolerance: Union[float, np.ndarray] = 1e-1, alpha: float = 0.05):

        self.num_cells = num_cells
        self.tolerance = tolerance
        self.alpha = alpha

        self.observed = dict()

    @classmethod
    def plot_sample_size(cls, num_cells: int, alpha: float = 0.05, **fit_kwargs) -> None:

        tolerance = np.exp(np.linspace(start=np.log(1e-2), stop=np.log(1e-1), num=100))

        estimator = cls(
            num_cells=num_cells,
            tolerance=tolerance,
            alpha=alpha
        ).fit(**fit_kwargs)

        confidence_level = (1 - estimator.alpha) * 100
        population_size = estimator.observed.get("population_size")

        required_sample_size = estimator.required_sample_size
        title_text = "Required sample size (cells: {} - {}conf: {}%)".format(
            num_cells,
            "" if population_size is None else "pop: {} - ".format(population_size),
            int(confidence_level)
        )

        plt.plot(tolerance, required_sample_size)
        plt.xlabel("Tolerance")
        plt.ylabel("Sample size")
        plt.title(title_text)
        plt.grid(True)

    @property
    def chi2(self) -> float:
        return scipy.stats.chi2.ppf(
            q=1 - self.alpha / self.num_cells,
            df=1
        )

    @property
    @requires_fit
    def required_sample_size(self) -> Union[int, np.ndarray]:

        population_size = self.observed.get("population_size")
        term = self.chi2 * self.proportions * (1. - self.proportions)

        denom = np.copy(term.reshape(-1, 1))

        tol = self.tolerance
        if not isinstance(self.tolerance, np.ndarray):
            tol = np.array([tol])

        nom = tol.reshape(1, -1) ** 2

        if population_size is not None:
            denom *= population_size
            nom = (population_size - 1) * nom + term.reshape(-1, 1)

        estimated_sizes = np.ceil((denom / nom).max(axis=0)).astype(int)
        if len(estimated_sizes) < 2:
            estimated_sizes = estimated_sizes.item()

        return estimated_sizes

    @property
    @requires_fit
    def sample_size(self) -> int:
        return self.observed["sample_size"]

    @property
    def proportions(self) -> np.ndarray:

        cell_frequencies = self.observed.get("cell_frequencies")
        sample_size = self.sample_size
        dtype = np.float32

        if cell_frequencies is not None:
            prop = np.array(cell_frequencies, dtype=dtype) / sample_size
        else:
            # No prior probability distribution available => use 'worst case'
            prop = np.array([0.5] + [0.5 / (self.num_cells - 1)] * (self.num_cells - 1), dtype=dtype)

        return prop

    @property
    def confidence_intervals(self):
        half_widths = self.get_confidence_half_widths()
        lower = self.proportions - half_widths
        upper = self.proportions + half_widths
        return np.clip(np.concatenate([
            lower.reshape(-1, 1), upper.reshape(-1, 1)
        ], axis=1), 0, 1)

    def get_confidence_half_widths(self) -> np.ndarray:

        half_widths = (
                self.chi2 * self.proportions * (1. - self.proportions) / self.sample_size
        )

        population_size = self.observed.get("population_size")
        if population_size is not None:
            half_widths *= (population_size - sample_size) / (population_size - 1)

        half_widths = np.sqrt(half_widths)

        return half_widths

    def fit(
            self,
            cell_frequencies: Optional[List[int]] = None,
            sample_size: Optional[int] = None,
            population_size: Optional[int] = None
    ) -> "MultinomialEstimator":

        _sample_size = sample_size

        if cell_frequencies is not None:

            _num_cells = len(cell_frequencies)
            if _num_cells != self.num_cells:
                raise ValueError("Number of cells of observation does not match the specified value!")

            _sample_size = sum(cell_frequencies)
            if sample_size is not None and sample_size != _sample_size:
                warnings.warn("You have specified a sample size that does not match your observation!")

            self.observed["cell_frequencies"] = cell_frequencies

        self.observed["sample_size"] = _sample_size

        if population_size is not None:
            self.observed["population_size"] = population_size

        return self
