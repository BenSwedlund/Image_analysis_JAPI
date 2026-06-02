from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import Mapping, Optional, Sequence

import h5py
import matplotlib.pyplot as plt
import pandas as pd

from analysis import process_image, analyze_cross_correlation
from viz import plot_cross_correlation_and_radial_profile


PathLike = str | Path
ClipLimits = tuple[float | None, float | None]


def _as_path(value: PathLike) -> Path:
    return value if isinstance(value, Path) else Path(value)


def _ensure_dir(path: Optional[PathLike]) -> Optional[Path]:
    if path is None:
        return None
    p = _as_path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _get_clip_limits(
    filename: str,
    clip_limits_by_name: Mapping[str, ClipLimits] | None,
) -> ClipLimits:
    if clip_limits_by_name is None:
        return (None, None)
    return clip_limits_by_name.get(filename, (None, None))


def main(
    ch1_files: Sequence[PathLike],
    ch2_files: Sequence[PathLike],
    analysis_output_dir: PathLike | None,
    sigma_um: float,
    max_corr_dist_um: float | int,
    pixel_size_um: float = 1 / 0.5507,
    nbins: int = 100,
    quantiles: Sequence[float] = (0.05, 0.95),
    conf_int_columns: Sequence[str] | None = ("Q_0.05", "Q_0.95"),
    save_all_data: bool = False,
    plot_output_dir: PathLike | None = None,
    clip_limits_by_name: Mapping[str, ClipLimits] | None = None,
    fmt: str = "png",
    dpi: int = 300,
) -> None:
    """
    Compute cross-correlation between paired images or channels and save results.
    """
    ch1_files = [_as_path(p) for p in ch1_files]
    ch2_files = [_as_path(p) for p in ch2_files]

    if len(ch1_files) != len(ch2_files):
        raise ValueError("The number of files in both channels must be the same.")

    out_dir = _ensure_dir(analysis_output_dir)
    plot_dir = _ensure_dir(plot_output_dir)

    saved_label = "all" if save_all_data else "summary"

    for ch1_f, ch2_f in zip(ch1_files, ch2_files):
        # Optional safety check: only pair matching stems
        if ch1_f.stem.replace("_mCherry", "") != ch2_f.stem.replace("_GFP", ""):
            print(f"Warning: stem mismatch for {ch1_f.name} and {ch2_f.name}. Skipping pair.")
            continue

        clip1 = _get_clip_limits(ch1_f.name, clip_limits_by_name)
        clip2 = _get_clip_limits(ch2_f.name, clip_limits_by_name)

        im1 = process_image(
            ch1_f,
            pixel_size_um=pixel_size_um,
            sigma_um=sigma_um,
            clip_limits=clip1,
            output_dir=None,
        )
        im2 = process_image(
            ch2_f,
            pixel_size_um=pixel_size_um,
            sigma_um=sigma_um,
            clip_limits=clip2,
            output_dir=None,
        )

        if im1.shape != im2.shape:
            print(
                f"Warning: image shapes do not match for {ch1_f.name} and {ch2_f.name}. Skipping pair."
            )
            continue

        results = analyze_cross_correlation(
            im1,
            im2,
            pixel_size_um=pixel_size_um,
            max_corr_dist_um=max_corr_dist_um,
            nbins=nbins,
            quantiles=quantiles,
            return_all_data=save_all_data,
        )

        if save_all_data:
            S, df, S_radial, D_radial = results
        else:
            S, df = results

        pair_name = f"{ch1_f.stem}_{ch2_f.stem}_cross_correlation"

        if out_dir is not None:
            out_file = out_dir / f"{pair_name}.h5"

            df.to_hdf(out_file, key="radial_profile", mode="w", format="fixed")
            with h5py.File(out_file, mode="a") as h:
                h.create_dataset("S", data=S)
                if save_all_data:
                    h.create_dataset("S_radial", data=S_radial)
                    h.create_dataset("D_radial", data=D_radial)

            print(f"Saved {saved_label} cross-correlation data to: {out_file.resolve()}")

        if plot_dir is not None:
            if conf_int_columns is None:
                conf_intervals = None
            else:
                missing = [c for c in conf_int_columns if c not in df.columns]
                if missing:
                    raise KeyError(f"Missing confidence interval columns: {missing}")
                conf_intervals = df.loc[:, list(conf_int_columns)].to_numpy().T

            plot_file = plot_dir / f"{date.today().isoformat()}_{pair_name}_plots.{fmt}"
            plot_cross_correlation_and_radial_profile(
                im1,
                im2,
                S=S,
                radial_distance_um=df["bin_centers_um"],
                S_mean=df["mean"],
                pixel_size_um=pixel_size_um,
                max_corr_dist_um=max_corr_dist_um,
                conf_intervals=conf_intervals,
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
        ch1_files=ch1_files,
        ch2_files=ch2_files,
        analysis_output_dir="/path/to/your/data/xcorr/analysis",
        sigma_um=0,
        max_corr_dist_um=1000,
        save_all_data=False,
        plot_output_dir="/path/to/your/data/xcorr/plots",
    )
