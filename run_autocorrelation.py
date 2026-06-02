from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import Optional, Sequence

import h5py
import matplotlib.pyplot as plt
import pandas as pd

from analysis import process_image, analyze_image_autocorrelation
from viz import plot_image_autocorrelation_and_radial_profile


PathLike = str | Path


def _as_path(value: PathLike) -> Path:
    return value if isinstance(value, Path) else Path(value)


def _ensure_dir(path: Optional[PathLike]) -> Optional[Path]:
    if path is None:
        return None
    p = _as_path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _validate_conf_int_columns(df: pd.DataFrame, columns: Sequence[str] | None) -> Optional[list[list[float]]]:
    if columns is None:
        return None

    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise KeyError(f"Missing confidence interval columns in dataframe: {missing}")

    return df.loc[:, list(columns)].to_numpy().T.tolist()


def main(
    image_files: Sequence[PathLike],
    analysis_output_dir: PathLike | None,
    pixel_size_um: float,
    sigma_um: float,
    max_corr_dist_um: float | int,
    nbins: int = 100,
    quantiles: Sequence[float] = (0.05, 0.95),
    conf_int_columns: Sequence[str] | None = ("Q_0.05", "Q_0.95"),
    save_all_data: bool = False,
    image_output_dir: PathLike | None = None,
    plot_output_dir: PathLike | None = None,
    clip_limits: Sequence[Sequence[float] | None] | None = None,
    fmt: str = "png",
    dpi: int = 300,
) -> None:
    """
    Compute image autocorrelation and radial profile for one or more images.

    Parameters
    ----------
    image_files
        Input image files to process.
    analysis_output_dir
        Directory for HDF5 output. If None, no analysis file is written.
    pixel_size_um
        Pixel size in microns.
    sigma_um
        Smoothing sigma passed to the image preprocessing step.
    max_corr_dist_um
        Maximum correlation distance in microns.
    nbins
        Number of radial bins.
    quantiles
        Quantiles used in the radial profile summary.
    conf_int_columns
        Columns to use for confidence intervals in the plot.
    save_all_data
        If True, save full autocorrelation arrays as well.
    image_output_dir
        Optional directory for processed images.
    plot_output_dir
        Optional directory for plots.
    clip_limits
        Optional per-image clip limits. Must match image_files length if provided.
    fmt
        Plot file format.
    dpi
        Plot resolution.
    """
    image_files = [_as_path(p) for p in image_files]

    if clip_limits is None:
        clip_limits = [None] * len(image_files)
    elif len(clip_limits) != len(image_files):
        raise ValueError("clip_limits and image_files must have the same length")

    out_dir = _ensure_dir(analysis_output_dir)
    image_dir = _ensure_dir(image_output_dir)
    plot_dir = _ensure_dir(plot_output_dir)

    saved_label = "all" if save_all_data else "summary"

    for image_file, clip_lims in zip(image_files, clip_limits):
        im = process_image(
            image_file,
            pixel_size_um=pixel_size_um,
            sigma_um=sigma_um,
            clip_limits=clip_lims,
            output_dir=image_dir,
        )

        results = analyze_image_autocorrelation(
            im,
            pixel_size_um=pixel_size_um,
            max_corr_dist_um=max_corr_dist_um,
            nbins=nbins,
            quantiles=quantiles,
            return_all_data=save_all_data,
        )

        if save_all_data:
            A, df, A_radial, D_radial = results
        else:
            A, df = results

        stem = image_file.stem

        if out_dir is not None:
            out_file = out_dir / f"{stem}_autocorrelation.h5"
            df.to_hdf(out_file, key="radial_profile", mode="w", format="fixed")

            with h5py.File(out_file, mode="a") as h:
                h.create_dataset("A", data=A)
                if save_all_data:
                    h.create_dataset("A_radial", data=A_radial)
                    h.create_dataset("D_radial", data=D_radial)

            print(f"Saved {saved_label} autocorrelation data to: {out_file.resolve()}")

        if plot_dir is not None:
            conf_intervals = _validate_conf_int_columns(df, conf_int_columns)

            plot_file = plot_dir / f"{date.today().isoformat()}_{stem}_radial_autocorrelation.{fmt}"
            plot_image_autocorrelation_and_radial_profile(
                im,
                A,
                radial_distance_um=df["bin_centers_um"],
                A_mean=df["mean"],
                pixel_size_um=pixel_size_um,
                max_corr_dist_um=max_corr_dist_um,
                conf_intervals=conf_intervals,
                title=stem,
                save=True,
                fname=plot_file,
                dpi=dpi,
            )
            plt.close()


if __name__ == "__main__":
    image_dir = Path("/path/to/your/data/")

    ch1_files = sorted(image_dir.glob("mCherry/*mCherry*.tif"))
    ch2_files = sorted(image_dir.glob("GFP/*GFP*.tif"))

    main(
        image_files=ch1_files + ch2_files,
        analysis_output_dir="/path/to/your/data/",
        pixel_size_um=1 / 0.5507,
        sigma_um=15,
        max_corr_dist_um=1000,
        nbins=100,
        save_all_data=True,
        # plot_output_dir="/path/to/your/data/plots",
    )
