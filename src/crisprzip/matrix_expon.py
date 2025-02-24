"""Perform fast and parallel matrix exponentiation."""

import numpy as np
from numba import njit
from scipy import linalg


@njit(cache=True)
def exponentiate_fast(matrix: np.ndarray, time: np.ndarray) -> np.ndarray:
    """Fast method to calculate exp(M * t), via diagonalizing the matrix M.
    Uses Numba's just-in-time compilation for performance optimization.

    Parameters
    ----------
    matrix : `numpy.ndarray`, (N, N)
        The input matrix to exponentiate.
    time : `numpy.ndarray`, (T,)
        A 1D array of time values at which to calculate the exponentiation.

    Returns
    -------
    exp_matrix : `numpy.ndarray`, (N, T, N)
        The resulting matrix exponentiation at each ``time`` step.
        Returns `None` if:
        - Matrix is singular
        - Matrix has complex eigenvalues
        - Overflow occurs in the exponential function
        - Resulting matrix contains negative terms.
    """

    # preallocate result matrix
    exp_matrix = np.zeros((matrix.shape[0], time.shape[0], matrix.shape[1]))

    # 1. diagonalize M = U D U_inv
    try:
        eigenvals, eigenvecs = np.linalg.eig(matrix)
        eigenvecs_inv = np.linalg.inv(eigenvecs)
    except Exception:
        return None  # Handles singular matrices
    if np.any(np.iscomplex(eigenvals)):
        return None  # Skip complex eigenvalue matrices

    for i in range(time.size):
        t = time[i]

        # 2. B = exp(Dt)
        exponent_matrix = eigenvals.real * t
        if np.any(exponent_matrix > 700.0):
            return None  # Avoids np.exp overflow
        b_matrix = np.diag(np.exp(exponent_matrix))

        # 3. exp(Mt) = U B U_inv = U exp(Dt) U_inv
        distr = eigenvecs @ b_matrix @ eigenvecs_inv
        if np.any(distr < 0.0):
            distr = np.maximum(distr, 0)
        exp_matrix[:, i, :] = distr

    return exp_matrix


@njit(cache=True)
def update_rate_matrix(ref_rate_matrix: np.ndarray,
                       on_rate: float) -> np.ndarray:
    """Update a reference rate matrix with a specific on-rate.

    Parameters
    ----------
    ref_rate_matrix : `numpy.ndarray`, (N, N)
        The reference transition rate matrix.
    on_rate : `float`
        The rate of transition to update in the matrix.

    Returns
    -------
    rate_matrix : `numpy.ndarray`, (N, N)
        A copy of the reference rate matrix with the ``on_rates`` value updated.
    """

    rate_matrix = ref_rate_matrix.copy()
    rate_matrix[0, 0] = -on_rate
    rate_matrix[1, 0] = on_rate
    return rate_matrix


@njit(cache=True)
def exponentiate_fast_var_onrate(ref_matrix: np.ndarray, time: float,
                                 on_rates: np.ndarray) -> np.ndarray:
    """Compute exp(M * t) for varying on-rates using diagonalization.

    Parameters
    ----------
    ref_matrix : `numpy.ndarray`, (N, N)
        The reference transition rate matrix.
    time : `float`
        The time at which the exponentiation is evaluated.
    on_rates : `numpy.ndarray`, (K,)
        A 1D array of varying on-rates to update the reference matrix.

    Returns
    -------
    exp_matrix : `numpy.ndarray`, (N, K, N)
        The resulting matrix exponentiation for each ``on_rates`` value.
        Returns `None` if diagonalization fails or produces invalid results.
    """

    exp_matrix = np.zeros(
        (ref_matrix.shape[0], on_rates.shape[0], ref_matrix.shape[1]))

    for i, k_on in enumerate(on_rates):
        rate_matrix = update_rate_matrix(ref_matrix, k_on)

        partial_result = exponentiate_fast(rate_matrix, np.array([time]))
        if partial_result is None:
            return None
        else:
            exp_matrix[:, i, :] = partial_result[:, 0, :]

    return exp_matrix


def exponentiate_iterative(matrix: np.ndarray, time: np.ndarray) -> np.ndarray:
    """Iteratively compute exp(M * t) using `scipy.linalg.expm`.

    Parameters
    ----------
    matrix : `numpy.ndarray`, (N, N)
        The input matrix to exponentiate.
    time : `numpy.ndarray`, (T,)
        A 1D array of time values at which to calculate the exponentiation.

    Returns
    -------
    exp_matrix : `numpy.ndarray`, (N, T, N)
        The resulting matrix exponentiation at each ``time`` step.
    """

    exp_matrix = np.zeros((matrix.shape[0], time.shape[0], matrix.shape[1]))
    for i in range(time.size):
        exp_matrix[:, i, :] = linalg.expm(matrix * time[i])
    return exp_matrix


def exponentiate_iterative_var_onrate(ref_matrix: np.ndarray, time: float,
                                      on_rates: np.ndarray) -> np.ndarray:
    """Iteratively compute exp(M * t) for varying on-rates using
    `scipy.linalg.expm`.

    Parameters
    ----------
    ref_matrix : `numpy.ndarray`, (N, N)
        The reference transition rate matrix.
    time : `float`
        The time at which the exponentiation is evaluated.
    on_rates : `numpy.ndarray`, (K,)
        A 1D array of varying on-rates to update the reference matrix.

    Returns
    -------
    exp_matrix : `numpy.ndarray`, (N, K, N)
        The resulting matrix exponentiation for each value in ``on-rates``.
    """

    exp_matrix = np.zeros(
        (ref_matrix.shape[0], on_rates.shape[0], ref_matrix.shape[1]))
    for i, k_on in enumerate(on_rates):
        rate_matrix = update_rate_matrix(ref_matrix, k_on)
        exp_matrix[:, i, :] = linalg.expm(rate_matrix * time)
    return exp_matrix
