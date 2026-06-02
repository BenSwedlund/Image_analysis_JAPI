from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from analysis import (
    manders_m1_m2,
    process_image,
    dice_sorenson_coeff,
    cohens_kappa,
)


PathLike = str | Path


def _as_path(value: PathLike) -> Path:
    return value if isinstance(value, Path) else Path(value)


def _ensure_dir(path: Optional[PathLike]) -> Optional[Path]:
    if path is None:
        return None
    p = _as_path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def plot_mask_overlap(
    mask1: np.ndarray,
    mask2: np.ndarray,
    pixel_size_um: float = 1 / 0.2753,
    figsize: tuple[int, int] = (10, 10),
    output_path: Optional[PathLike] = None,
    fmt: str = "png",
    dpi: int = 300,
) -> None:
    """
    Plot two binary masks in red/green with overlap in yellow.
    """
    import matplotlib.pyplot as plt

    overlap = np.zeros((*mask1.shape, 3), dtype=np.uint8)
    mask1_pos = mask1 > 0
    mask2_pos = mask2 > 0

    overlap[mask1_pos] = [255, 0, 0]
    overlap[mask2_pos] = [0, 255, 0]
    overlap[mask1_pos & mask2_pos] = [255, 255, 0]

    plt.figure(figsize=figsize)
    plt.imshow(overlap)
    plt.axis("equal")
    plt.xlabel("Distance (um)")
    plt.ylabel("Distance (um)")

    step = max(1, int(10 / pixel_size_um))
    x_ticks = np.arange(0, mask1.shape[1], step)
    y_ticks = np.arange(0, mask1.shape[0], step)
    plt.xticks(x_ticks, (x_ticks * pixel_size_um).astype(int))
    plt.yticks(y_ticks, (y_ticks * pixel_size_um).astype(int))

    if output_path is not None:
        out = _as_path(output_path)
        print(f"Saving overlap plot to: {out.resolve()}")
        plt.savefig(out, format=fmt, dpi=dpi, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def analyze_overlap(
    ch1_file: PathLike,
    ch2_file: PathLike,
    pixel_size_um: float = 1 / 0.2753,
) -> tuple[np.ndarray, np.ndarray, float, float, float, float, float, float, float, float]:
    """
    Compute overlap metrics for two binary masks.

    Returns
    -------
    mask1, mask2, kappa, m1, m2, dice, mask1_mean, mask2_mean, overlap_mean, expected_overlap
    """
    ch1_file = _as_path(ch1_file)
    ch2_file = _as_path(ch2_file)

    mask1 = process_image(
        ch1_file,
        pixel_size_um=pixel_size_um,
        sigma_um=0,
        clip_limits=None,
        output_dir=None,
    )
    mask2 = process_image(
        ch2_file,
        pixel_size_um=pixel_size_um,
        sigma_um=0,
        clip_limits=None,
        output_dir=None,
    )

    m1, m2 = manders_m1_m2(mask1, mask2)
    dice = dice_sorenson_coeff(mask1, mask2)
    kappa = cohens_kappa(mask1, mask2)

    mask1_pos = mask1 > 0
    mask2_pos = mask2 > 0

    mask1_mean = mask1_pos.mean()
    mask2_mean = mask2_pos.mean()
    overlap_mean = (mask1_pos & mask2_pos).mean()

    # Expected overlap if the masks were independent
    expected_overlap = mask1_mean * mask2_mean

    return (
        mask1,
        mask2,
        kappa,
        m1,
        m2,
        dice,
        mask1_mean,
        mask2_mean,
        overlap_mean,
        expected_overlap,
    )


def main(
    ch1_files: Sequence[PathLike],
    ch2_files: Sequence[PathLike],
    analysis_output_dir: PathLike | None,
    pixel_size_um: float = 1 / 0.2753,
    plot_output_dir: PathLike | None = None,
    fmt: str = "png",
    dpi: int = 300,
) -> None:
    ch1_files = [_as_path(p) for p in ch1_files]
    ch2_files = [_as_path(p) for p in ch2_files]

    if len(ch1_files) != len(ch2_files):
        raise ValueError("The number of files in both channels must be the same.")

    out_dir = _ensure_dir(analysis_output_dir)
    plot_dir = _ensure_dir(plot_output_dir)

    n_samples = len(ch1_files)

    mask1_means = np.zeros(n_samples)
    mask2_means = np.zeros(n_samples)
    overlap_means = np.zeros(n_samples)

    m1s = np.zeros(n_samples)
    m2s = np.zeros(n_samples)
    kappas = np.zeros(n_samples)
    dices = np.zeros(n_samples)

    expected_overlaps = np.zeros(n_samples)

    for i, (ch1_file, ch2_file) in enumerate(zip(ch1_files, ch2_files)):
        (
            im1,
            im2,
            kappa,
            m1,
            m2,
            dice,
            mask1_mean,
            mask2_mean,
            overlap_mean,
            expected_overlap,
        ) = analyze_overlap(ch1_file, ch2_file, pixel_size_um)

        mask1_means[i] = mask1_mean
        mask2_means[i] = mask2_mean
        overlap_means[i] = overlap_mean

        kappas[i] = kappa
        m1s[i] = m1
        m2s[i] = m2
        dices[i] = dice
        expected_overlaps[i] = expected_overlap

        print(f"Overlap analysis for {ch1_file.name} and {ch2_file.name}:")
        print(f" -- mask1_mean: {mask1_mean:.3f}")
        print(f" -- mask2_mean: {mask2_mean:.3f}")
        print(f" -- overlap_mean: {overlap_mean:.3f}")
        print(f" --  M1: {m1:.3f}")
        print(f" --  M2: {m2:.3f}")
        print(f" -- D-S: {dice:.3f}")
        print(f" -- expected overlap: {expected_overlap:.3f}")

        if plot_dir is not None:
            today = date.today().isoformat()
            plot_path = plot_dir / f"{today}_{ch1_file.stem}_{ch2_file.stem}_overlap.{fmt}"
            plot_mask_overlap(
                im1,
                im2,
                pixel_size_um=pixel_size_um,
                output_path=plot_path,
                fmt=fmt,
                dpi=dpi,
            )

        print()

    m1s_excess = m1s - mask2_means
    m2s_excess = m2s - mask1_means
    # Excess overlap above random expectation
    dices_excess = dices - (
        2 * mask1_means * mask2_means / (mask1_means + mask2_means)
    )

    # Colocalization coefficient (observed / expected Dice)
    expected_dices = (
        2 * mask1_means * mask2_means / (mask1_means + mask2_means)
    )

    dices_excess_ratio = np.divide(
        dices,
        expected_dices,
        out=np.full_like(dices, np.nan),
        where=expected_dices > 0,
    )

    if out_dir is not None:
        df = pd.DataFrame(
            {
                "ch1_file": [f.name for f in ch1_files],
                "ch2_file": [f.name for f in ch2_files],
                "mask1_mean": mask1_means,
                "mask2_mean": mask2_means,
                "overlap_mean": overlap_means,
                "expected_overlap": expected_overlaps,
                "Cohen_kappa": kappas,
                "M1_coeff": m1s,
                "M2_coeff": m2s,
                "Dice_coeff": dices,
                "Excess_overlap": dices_excess,
                "Colocalization_coeff": dices_excess_ratio,
                "M1_excess": m1s_excess,
                "M2_excess": m2s_excess,
            }
        )

        today = date.today().isoformat()
        output_path = out_dir / f"{today}_overlap_analysis_results.csv"
        print(f"Saving overlap analysis results to: {output_path.resolve()}")
        df.to_csv(output_path, index=False)


if __name__ == "__main__":
    image_dir = Path("/path/to/your/data/")
    ch1_files = sorted(image_dir.glob("*mCherry*.tif"))
    ch2_files = sorted(image_dir.glob("*GFP*.tif"))

    main(
        ch1_files=ch1_files,
        ch2_files=ch2_files,
        analysis_output_dir="/path/to/your/data/",
    )
