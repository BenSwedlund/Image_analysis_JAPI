from __future__ import annotations

from math import floor
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import fft

from image import blur_image, load_image


__all__ = [
    "manders_m1_m2",
    "cohens_kappa",
    "dice_sorenson_coeff",
    "sort_2d_array_by_distance",
    "radial_profile",
    "autocorrelation_2d",
    "cross_correlation_2d",
    "process_image",
    "analyze_radial_profile",
    "analyze_image_autocorrelation",
    "analyze_cross_correlation",
]


PathLike = str | Path


def _as_path(value: PathLike) -> Path:
    return value if isinstance(value, Path) else Path(value)


def _ensure_dir(path: PathLike | None) -> Path | None:
    if path is None:
        return None
    p = _as_path(path)
    if not p.exists():
        raise ValueError(f"Directory does not exist: {p}")
    if not p.is_dir():
        raise ValueError(f"Path is not a directory: {p}")
    return p


def manders_m1_m2(x1: np.ndarray, x2: np.ndarray) -> tuple[float, float]:
    """
    Compute Manders M1 and M2 for non-negative arrays.

    Here, colocalization is defined on pixels where both arrays are > 0.
    """
    if x1.shape != x2.shape:
        raise ValueError("Input images must have the same shape.")
    if np.any(x1 < 0) or np.any(x2 < 0):
        raise ValueError("Input images must have non-negative values.")

    x1_mask = x1 > 0
    x2_mask = x2 > 0
    coloc = x1_mask & x2_mask

    x1_total = x1.sum()
    x2_total = x2.sum()
    if x1_total == 0 or x2_total == 0:
        return np.nan, np.nan

    m1 = x1[coloc].sum() / x1_total
    m2 = x2[coloc].sum() / x2_total
    return float(m1), float(m2)


def cohens_kappa(x1: np.ndarray, x2: np.ndarray) -> float:
    """Compute Cohen's kappa for two binary masks."""
    if x1.shape != x2.shape:
        raise ValueError("Input images must have the same shape.")

    x1_mask = x1 > 0
    x2_mask = x2 > 0

    if np.array_equal(x1_mask, x2_mask):
        return 1.0

    n = x1_mask.size

    a = (x1_mask & x2_mask).sum()
    b = x1_mask.sum() - a
    c = x2_mask.sum() - a
    d = n - (a + b + c)

    p0 = (a + d) / n
    pe = ((a + b) * (a + c) + (c + d) * (b + d)) / n**2
    if pe == 1:
        return np.nan

    return float((p0 - pe) / (1 - pe))


def dice_sorenson_coeff(x1: np.ndarray, x2: np.ndarray) -> float:
    """Compute Dice-Sørensen coefficient for two non-negative arrays."""
    if x1.shape != x2.shape:
        raise ValueError("Input images must have the same shape.")
    if np.any(x1 < 0) or np.any(x2 < 0):
        raise ValueError("Input images must have non-negative values.")

    x1_mask = x1 > 0
    x2_mask = x2 > 0
    denom = x1_mask.sum() + x2_mask.sum()
    if denom == 0:
        return np.nan

    intersection = (x1_mask & x2_mask).sum()
    return float(2 * intersection / denom)


def sort_2d_array_by_distance(
    X: np.ndarray,
    center: tuple[float, float] | None = None,
    scale: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Flatten a 2D array and sort values by distance from a center point.
    """
    if X.ndim != 2:
        raise ValueError("X must be 2-dimensional.")

    m, n = X.shape
    yy, xx = np.mgrid[:m, :n]
    xx = xx.ravel()
    yy = yy.ravel()

    if center is None:
        y_c = m / 2
        x_c = n / 2
    else:
        y_c, x_c = center

    distances = np.sqrt((xx - x_c) ** 2 + (yy - y_c) ** 2)
    if scale is not None:
        distances = distances * scale

    order = np.argsort(distances)
    return X.ravel()[order], distances[order]


def radial_profile(
    X: np.ndarray,
    nbins: int | None = None,
    bin_edges_px: np.ndarray | None = None,
    center: tuple[float, float] | None = None,
    max_dist_px: float | None = None,
    mean: bool = True,
    quantiles: Sequence[float] = (0.5,),
    normalize: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute radial summary statistics of a 2D image.
    Returns bin edges, number of observations per bin, and metrics.
    """
    if not mean and len(quantiles) == 0:
        raise ValueError("Pass at least one quantile or set mean=True.")

    if X.ndim != 2:
        raise ValueError("X must be 2-dimensional.")

    m, n = X.shape
    if max_dist_px is None:
        max_dist_px = min(m, n) / 2

    if bin_edges_px is None:
        if nbins is None:
            nbins = max(1, floor(np.sqrt(X.size)))
        bin_edges_px = np.linspace(0, max_dist_px, nbins + 1)
    else:
        nbins = len(bin_edges_px) - 1

    X_flat, D_flat_px = sort_2d_array_by_distance(X, center=center, scale=None)
    keep = D_flat_px <= max_dist_px
    X_flat = X_flat[keep]
    D_flat_px = D_flat_px[keep]

    if normalize:
        max_abs = np.abs(X_flat).max()
        if max_abs != 0:
            X_flat = X_flat / max_abs

    bin_indices = np.digitize(D_flat_px, bins=bin_edges_px) - 1

    n_metrics = len(quantiles) + int(mean)
    metrics = np.full((nbins, n_metrics), np.nan, dtype=float)
    n_obs = np.zeros(nbins, dtype=int)

    for i in range(nbins):
        mask = bin_indices == i
        n_obs[i] = int(mask.sum())
        if n_obs[i] == 0:
            continue

        x_bin = X_flat[mask]
        q_vals = np.quantile(x_bin, q=quantiles) if len(quantiles) > 0 else np.array([])
        if len(quantiles) > 0:
            metrics[i, : len(quantiles)] = q_vals
        if mean:
            metrics[i, -1] = x_bin.mean()

    return bin_edges_px, n_obs, metrics


def autocorrelation_2d(X: np.ndarray) -> np.ndarray:
    """Compute the normalized 2D autocorrelation of an image."""
    X = np.asarray(X)
    F = fft.fft2(X - X.mean())
    P = np.abs(F) ** 2
    acf = fft.ifft2(P).real
    acf = fft.fftshift(acf)

    max_abs = np.abs(acf).max()
    if max_abs != 0:
        acf = acf / max_abs

    return acf


def cross_correlation_2d(im1: np.ndarray, im2: np.ndarray) -> np.ndarray:
    """Compute the normalized 2D cross-correlation of two images."""
    if im1.shape != im2.shape:
        raise ValueError("Input images must have the same shape.")

    im1 = np.asarray(im1, dtype=float)
    im2 = np.asarray(im2, dtype=float)

    im1_std = im1.std()
    im2_std = im2.std()
    if im1_std == 0 or im2_std == 0:
        raise ValueError("Input images must have non-zero standard deviation.")

    im1 = (im1 - im1.mean()) / im1_std
    im2 = (im2 - im2.mean()) / im2_std

    F1 = fft.fft2(im1)
    F2 = fft.fft2(im2)
    xcf = fft.ifft2(F1 * np.conj(F2)).real
    xcf = fft.fftshift(xcf)

    max_abs = np.abs(xcf).max()
    if max_abs != 0:
        xcf = xcf / max_abs

    return xcf


def process_image(
    file: PathLike,
    pixel_size_um: float,
    sigma_um: float | None,
    clip_limits: tuple[float | None, float | None] | None = None,
    output_dir: PathLike | None = None,
) -> np.ndarray:
    """
    Load and optionally blur an image. Optionally save the processed image.
    """
    out_dir = _ensure_dir(output_dir)
    imfile = _as_path(file)

    im = load_image(str(imfile), as_gray=True)
    if sigma_um is not None and sigma_um > 0:
        im = blur_image(
            im,
            sigma_um=sigma_um,
            pixel_size_um=pixel_size_um,
            clip=clip_limits,
        )

    if out_dir is not None:
        from skimage.io import imsave

        out_file = out_dir / f"{imfile.stem}_processed.tif"
        imsave(str(out_file), im)
        print(f"Saved processed image to: {out_file.resolve()}")

    return im


def analyze_radial_profile(
    Z: np.ndarray,
    max_corr_dist_um: float | int,
    pixel_size_um: float,
    nbins: int = 100,
    quantiles: Sequence[float] = (0.1, 0.5, 0.9),
    mean: bool = True,
) -> pd.DataFrame:
    """
    Compute radial summary statistics and return a tidy DataFrame.
    """
    bin_edges_px, n_obs, metrics = radial_profile(
        Z,
        nbins=nbins,
        max_dist_px=max_corr_dist_um / pixel_size_um,
        quantiles=quantiles,
        mean=mean,
    )

    bin_edges_um = bin_edges_px * pixel_size_um
    bin_centers_um = bin_edges_um[:-1] + np.diff(bin_edges_um) / 2

    col_names = [f"Q_{q:.2f}" for q in quantiles]
    if mean:
        col_names.append("mean")

    df = pd.DataFrame(
        {
            "bin_centers_um": bin_centers_um,
            "n_observations": n_obs,
            **{name: metrics[:, i] for i, name in enumerate(col_names)},
        }
    )
    return df


def analyze_image_autocorrelation(
    im: np.ndarray,
    pixel_size_um: float,
    max_corr_dist_um: float | int,
    nbins: int = 100,
    quantiles: Sequence[float] = (0.1, 0.5, 0.9),
    return_all_data: bool = False,
) -> tuple:
    """
    Compute autocorrelation and its radial summary statistics.
    """
    A = autocorrelation_2d(im)
    df = analyze_radial_profile(
        A,
        max_corr_dist_um=max_corr_dist_um,
        pixel_size_um=pixel_size_um,
        nbins=nbins,
        quantiles=quantiles,
        mean=True,
    )

    if return_all_data:
        return A, df, *sort_2d_array_by_distance(A, scale=pixel_size_um)
    return A, df


def analyze_cross_correlation(
    im1: np.ndarray,
    im2: np.ndarray,
    pixel_size_um: float,
    max_corr_dist_um: float | int,
    nbins: int = 100,
    quantiles: Sequence[float] = (0.1, 0.5, 0.9),
    return_all_data: bool = False,
) -> tuple:
    """
    Compute cross-correlation and its radial summary statistics.
    """
    S = cross_correlation_2d(im1, im2)
    df = analyze_radial_profile(
        S,
        max_corr_dist_um=max_corr_dist_um,
        pixel_size_um=pixel_size_um,
        nbins=nbins,
        quantiles=quantiles,
        mean=True,
    )

    if return_all_data:
        return S, df, *sort_2d_array_by_distance(S, scale=pixel_size_um)
    return S, df
