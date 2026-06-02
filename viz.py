from __future__ import annotations

from math import ceil
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

from analysis import sort_2d_array_by_distance

__all__ = [
    "plot_image_autocorrelation_and_radial_profile",
    "plot_cross_correlation_and_radial_profile",
]


def _hide_top_right_spines(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _validate_2d_arrays(x1: np.ndarray, x2: np.ndarray) -> None:
    if x1.ndim != 2 or x2.ndim != 2:
        raise ValueError("Both input arrays must be 2-dimensional.")
    if x1.shape != x2.shape:
        raise ValueError("Input arrays must have matching shapes.")


def _save_or_show(
    fig: plt.Figure,
    save: bool,
    fname: Optional[Path | str],
    fmt: str,
    dpi: int,
) -> None:
    if save:
        if fname is None:
            raise ValueError("fname must be provided when save=True.")
        out = Path(fname)
        fig.savefig(out, format=fmt, dpi=dpi, bbox_inches="tight")
        print(f"Plot saved to: {out.resolve()}")
    else:
        plt.show()


def plot_image_autocorrelation_and_radial_profile(
    X: np.ndarray,
    A: np.ndarray,
    radial_distance_um: Sequence[float] | np.ndarray,
    A_mean: Sequence[float] | np.ndarray,
    pixel_size_um: float,
    max_corr_dist_um: float,
    conf_intervals: Optional[np.ndarray] = None,
    title: str = "Preprocessed Image",
    fig: Optional[plt.Figure] = None,
    save: bool = False,
    fname: Optional[Path | str] = None,
    fmt: str = "png",
    dpi: int = 300,
) -> None:
    """
    Plot an image, its autocorrelation, and the radial autocorrelation profile.
    """
    if fig is None:
        fig = plt.figure(figsize=(15, 5))

    X = np.asarray(X)
    A = np.asarray(A)
    radial_distance_um = np.asarray(radial_distance_um)
    A_mean = np.asarray(A_mean)

    m, n = X.shape
    height_um = m * pixel_size_um
    width_um = n * pixel_size_um

    w = ceil(max_corr_dist_um / pixel_size_um)
    m_center = m // 2
    n_center = n // 2
    A_trunc = A[m_center - w : m_center + w, n_center - w : n_center + w]

    ax = fig.add_subplot(1, 3, 1)
    ax.imshow(X, cmap=plt.cm.copper, interpolation="bilinear", extent=(0, width_um, 0, height_um))
    ax.set_title(title)
    ax.set_xlabel("Distance (um)")
    ax.set_ylabel("Distance (um)")
    ax.axis("equal")
    ax.set_xticks(np.linspace(0, width_um, 5))
    ax.set_xticklabels(np.round(np.linspace(0, width_um, 5), 1))
    ax.set_yticks(np.linspace(0, height_um, 5))
    ax.set_yticklabels(np.round(np.linspace(0, height_um, 5), 1))
    _hide_top_right_spines(ax)

    ax = fig.add_subplot(1, 3, 2)
    im = ax.imshow(
        A_trunc,
        cmap=plt.cm.RdBu_r,
        vmin=-1,
        vmax=1,
        extent=(-max_corr_dist_um, max_corr_dist_um, -max_corr_dist_um, max_corr_dist_um),
    )
    ax.set_title("Autocorrelation")
    ax.set_xlabel("X displacement (um)")
    ax.set_ylabel("Y displacement (um)")
    ax.set_xticks(np.linspace(-max_corr_dist_um, max_corr_dist_um, 5))
    ax.set_xticklabels(np.round(np.linspace(-max_corr_dist_um, max_corr_dist_um, 5), 1))
    ax.set_yticks(np.linspace(-max_corr_dist_um, max_corr_dist_um, 5))
    ax.set_yticklabels(np.round(np.linspace(-max_corr_dist_um, max_corr_dist_um, 5), 1))
    fig.colorbar(im, ax=ax, shrink=0.7)
    _hide_top_right_spines(ax)

    ax = fig.add_subplot(1, 3, 3)
    ax.plot(radial_distance_um, A_mean, lw=1)

    min_yval = min(float(np.min(A_mean)), 0.0)
    if conf_intervals is not None:
        conf_intervals = np.asarray(conf_intervals)
        ax.fill_between(
            radial_distance_um,
            conf_intervals[0],
            conf_intervals[1],
            alpha=0.3,
            color="gray",
        )
        min_yval = min(min_yval, float(np.min(conf_intervals)))

    ax.set_title("Average radial autocorrelation")
    ax.set_xlabel("Distance (um)")
    ax.set_ylabel("Avg. autocorrelation")
    ax.set_ylim(min_yval, None)

    fig.tight_layout()
    _save_or_show(fig, save, fname, fmt, dpi)


def plot_cross_correlation_and_radial_profile(
    X1: np.ndarray,
    X2: np.ndarray,
    S: np.ndarray,
    radial_distance_um: Sequence[float] | np.ndarray,
    S_mean: Sequence[float] | np.ndarray,
    pixel_size_um: float,
    max_corr_dist_um: float,
    fig: Optional[plt.Figure] = None,
    conf_intervals: Optional[np.ndarray] = None,
    save: bool = False,
    fname: Optional[Path | str] = None,
    fmt: str = "png",
    dpi: int = 300,
    figsize: tuple[int, int] = (10, 10),
) -> None:
    """
    Plot two images, their cross-correlation, and the radial cross-correlation profile.
    """
    _validate_2d_arrays(X1, X2)

    if fig is None:
        fig = plt.figure(figsize=figsize)

    X1 = np.asarray(X1)
    X2 = np.asarray(X2)
    S = np.asarray(S)
    radial_distance_um = np.asarray(radial_distance_um)
    S_mean = np.asarray(S_mean)

    m, n = X1.shape
    height_um = m * pixel_size_um
    width_um = n * pixel_size_um
    height_mm = height_um / 1000.0
    width_mm = width_um / 1000.0

    w = ceil(max_corr_dist_um / pixel_size_um)
    m_center = m // 2
    n_center = n // 2
    S_trunc = S[m_center - w : m_center + w, n_center - w : n_center + w]

    S_trunc_flat, D_flat_px = sort_2d_array_by_distance(S_trunc)
    norm_region = np.abs(S_trunc_flat[D_flat_px <= w])
    norm_factor = norm_region.max() if norm_region.size else 1.0
    if norm_factor == 0:
        norm_factor = 1.0
    S_trunc = S_trunc / norm_factor

    ax = fig.add_subplot(2, 2, 1)
    x1_std = X1.std()
    x2_std = X2.std()
    x1_rescaled = (X1 - X1.mean()) / x1_std if x1_std != 0 else X1 - X1.mean()
    x2_rescaled = (X2 - X2.mean()) / x2_std if x2_std != 0 else X2 - X2.mean()

    xyprod = x1_rescaled * x2_rescaled

    x1_min = X1.min()
    x1_max = X1.max()
    x2_min = X2.min()
    x2_max = X2.max()

    x1_norm = (X1 - x1_min) / (x1_max - x1_min) if x1_max != x1_min else np.zeros_like(X1)
    x2_norm = (X2 - x2_min) / (x2_max - x2_min) if x2_max != x2_min else np.zeros_like(X2)
    blue = np.clip(xyprod, 0, None)
    blue = blue / blue.max() if blue.max() != 0 else blue

    ax.imshow(np.stack([x1_norm, x2_norm, blue], axis=-1), extent=(0, width_mm, 0, height_mm))
    ax.set_title("Channel 1 (R), Channel 2 (G), positive correlation (B)")
    ax.set_xlabel("Distance (mm)")
    ax.set_ylabel("Distance (mm)")
    ax.axis("equal")
    ax.set_xticks(np.linspace(0, width_mm, 5))
    ax.set_xticklabels(np.round(np.linspace(0, width_mm, 5), 1))
    ax.set_yticks(np.linspace(0, height_mm, 5))
    ax.set_yticklabels(np.round(np.linspace(0, height_mm, 5), 1))
    _hide_top_right_spines(ax)

    ax = fig.add_subplot(2, 2, 2)
    xyprod = xyprod / np.abs(xyprod).max() if np.abs(xyprod).max() != 0 else xyprod
    im = ax.imshow(
        xyprod,
        cmap=plt.cm.RdBu_r,
        vmin=-1,
        vmax=1,
        extent=(0, width_mm, 0, height_mm),
    )
    fig.colorbar(im, ax=ax, shrink=0.7)
    ax.set_title("Product of mean-centered channels")
    ax.set_xlabel("Distance (mm)")
    ax.set_ylabel("Distance (mm)")
    _hide_top_right_spines(ax)

    ax = fig.add_subplot(2, 2, 3)
    im = ax.imshow(
        S_trunc,
        cmap=plt.cm.RdBu_r,
        vmin=-1,
        vmax=1,
        extent=(-max_corr_dist_um, max_corr_dist_um, -max_corr_dist_um, max_corr_dist_um),
    )
    ax.set_title("Cross-correlation")
    ax.set_xlabel("X displacement (um)")
    ax.set_ylabel("Y displacement (um)")
    ax.set_xticks(np.linspace(-max_corr_dist_um, max_corr_dist_um, 5))
    ax.set_xticklabels(np.round(np.linspace(-max_corr_dist_um, max_corr_dist_um, 5), 1))
    ax.set_yticks(np.linspace(-max_corr_dist_um, max_corr_dist_um, 5))
    ax.set_yticklabels(np.round(np.linspace(-max_corr_dist_um, max_corr_dist_um, 5), 1))
    fig.colorbar(im, ax=ax, shrink=0.7)
    _hide_top_right_spines(ax)

    ax = fig.add_subplot(2, 2, 4)
    ax.plot(radial_distance_um, S_mean, lw=1)

    min_yval = min(float(np.min(S_mean)), 0.0)
    if conf_intervals is not None:
        conf_intervals = np.asarray(conf_intervals)
        ax.fill_between(
            radial_distance_um,
            conf_intervals[0],
            conf_intervals[1],
            alpha=0.3,
            color="gray",
        )
        min_yval = min(min_yval, float(np.min(conf_intervals)))

    ax.set_xlabel("Distance (um)")
    ax.set_ylabel("Mean cross-correlation")
    ax.set_ylim(min_yval, None)

    fig.tight_layout()
    _save_or_show(fig, save, fname, fmt, dpi)
