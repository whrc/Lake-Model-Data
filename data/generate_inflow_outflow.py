#!/usr/bin/env python3
"""
Build LAKE-model tributary inflow / outflow ASCII files from a meteo CSV and a driver file.

- Meteo file: same format as LAKE_NSF/meteo/*.dat (Year, Month, Day, Temp in K, Precip in m/s, ...).
- Driver file: LAKE driver; initial liquid water depth is read from keyword ``h10`` (m).
- Inflow temperature follows air temperature (meteo Temp converted to °C).
- Other inflow fields are static scalars (optional arguments).
- Outflow ``level`` (m): outside the warm season, fixed at ``h10``. Each calendar year, the
  leaky bucket runs only on the **longest consecutive stretch** with air temperature above a
  threshold (default 5 °C), starting from ``h10`` again each year.   Modeled depth is capped at
  ``h10`` (never above initial level), then ``level`` is a 5-day centered rolling mean.
  Hamon PET as in leaky_bucket.ipynb.

Use ``--plot`` to write a PNG next to this script (``<meteo_stem>_leaky_bucket_timeseries.png``).

Column order matches gen_in_outflows_file.ipynb (no CSV headers):
  Inflow:  Date, U, depth, width, temp, sal, Ux, Uy, DOC, POC, DIC, CH4, O2
  Outflow: Date, U, depth0, depth, width, level
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# Plots (saved next to this script; no GUI required)
# ---------------------------------------------------------------------------


def _safe_plot_stem(name: str) -> str:
    import re

    n = re.sub(r"[^\w\-.]+", "_", (name or "plot").strip())
    n = n.strip("_")[:100]
    return n or "leaky_bucket"


def save_leaky_bucket_timeseries_plot(
    dates: pd.Series,
    precip_m_per_day: np.ndarray,
    temp_c: np.ndarray,
    level_m: np.ndarray,
    *,
    title: str = "",
    output_path: Optional[Union[str, Path]] = None,
) -> Path:
    """
    Save a ``leaky_bucket.ipynb``-style figure (precip + water level) as PNG.

    Default path: ``<this_script_dir>/<title>_leaky_bucket_timeseries.png``
    Uses the Agg backend (works without a display).
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Patch

    dates = pd.to_datetime(dates)
    plot_df = pd.DataFrame(
        {
            "Date": dates,
            "Precip_mm": np.asarray(precip_m_per_day, dtype=float) * 1000.0,
            "Temp": np.asarray(temp_c, dtype=float),
            "depth_m": np.asarray(level_m, dtype=float),
        }
    )

    colors = np.where(plot_df["Temp"].to_numpy(dtype=float) > 0, "darkblue", "lightblue")
    rain_patch = Patch(color="darkblue", label="Rain (T > 0°C)", alpha=0.5)
    snow_patch = Patch(color="lightblue", label="Snow (T ≤ 0°C)", alpha=0.5)

    site_title = title or "Lake"
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.5))
    ax.set_title(f"Leaky Bucket Model: {site_title}")

    sns.lineplot(
        data=plot_df,
        x="Date",
        y="depth_m",
        ax=ax,
        color="red",
        label="Modeled water level",
        linewidth=0.8,
    )
    ax.set_ylabel("Water level (m)")
    ax.set_xlabel("Date")

    ax_precip = ax.twinx()
    ax_precip.bar(
        plot_df["Date"],
        plot_df["Precip_mm"],
        width=1.0,
        color=colors,
        linewidth=0.0,
        align="center",
        alpha=0.5,
    )
    ax_precip.set_ylabel("Precipitation (mm/day)")
    ax_precip.set_xlim(ax.get_xlim())

    h, lab = ax.get_legend_handles_labels()
    h.extend([rain_patch, snow_patch])
    lab.extend(["Rain (T > 0°C)", "Snow (T ≤ 0°C)"])
    ax.legend(h, lab, frameon=False, loc="best")

    fig.autofmt_xdate()
    plt.tight_layout()

    if output_path is None:
        out = _SCRIPT_DIR / f"{_safe_plot_stem(site_title)}_leaky_bucket_timeseries.png"
    else:
        out = Path(output_path)

    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out.resolve()


# ---------------------------------------------------------------------------
# Driver parsing (Situkuyuk_driver.dat style)
# ---------------------------------------------------------------------------


def _read_driver_keyword(driver_path: Union[str, Path], keyword: str) -> float:
    """Return the numeric value after ``keyword`` on a non-comment line."""
    path = Path(driver_path)
    with path.open(encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[0] == keyword:
                return float(parts[1])
    raise ValueError(f"Keyword '{keyword}' not found in {path}")


def read_initial_lake_depth_m(driver_path: Union[str, Path]) -> float:
    """Initial liquid water thickness from driver ``h10`` (m)."""
    return _read_driver_keyword(driver_path, "h10")


def read_latitude_deg(driver_path: Union[str, Path]) -> float:
    """Latitude from driver ``phi`` (deg N), for Hamon PET."""
    return _read_driver_keyword(driver_path, "phi")


# ---------------------------------------------------------------------------
# Leaky bucket + PET (from leaky_bucket.ipynb; PET fixed to use explicit dates)
# ---------------------------------------------------------------------------


def pond_water_balance(
    params: Tuple[float, float, float],
    precip_m: np.ndarray,
    pet_m: np.ndarray,
    initial_depth_m: float,
) -> np.ndarray:
    """
    Daily depth (m). params = (catchment_runoff_ratio, seepage_m_day, et_scale).
    inflow = precip * (1 + c_ratio); outflow = seepage + PET * et_scale.
    Depth is clamped to [0, initial_depth_m] each day (cannot exceed initial level).
    """
    c_ratio, seepage_m_day, et_scale = params
    h_max = float(initial_depth_m)
    depths = np.zeros(len(precip_m), dtype=float)
    current = h_max
    for t in range(len(precip_m)):
        inflow = precip_m[t] * (1.0 + c_ratio)
        outflow = seepage_m_day + pet_m[t] * et_scale
        current = max(0.0, min(h_max, current + inflow - outflow))
        depths[t] = current
    return depths


def longest_consecutive_run_indices(mask: np.ndarray) -> Optional[Tuple[int, int]]:
    """
    Return inclusive (start, end) row indices for the longest contiguous run of True values.

    If ``mask`` is empty or has no True values, return None. On ties, the earliest run wins.
    """
    mask = np.asarray(mask, dtype=bool)
    n = mask.size
    best: Optional[Tuple[int, int]] = None
    best_len = 0
    i = 0
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while j < n and mask[j]:
            j += 1
        length = j - i
        if length > best_len:
            best_len = length
            best = (i, j - 1)
        i = j
    return best


def level_from_seasonal_leaky_bucket(
    meteo_df: pd.DataFrame,
    precip_m: np.ndarray,
    pet_m: np.ndarray,
    leaky_params: Tuple[float, float, float],
    initial_depth_m: float,
    *,
    warm_temp_threshold_c: float = 5.0,
) -> np.ndarray:
    """
    Per calendar year: find the longest consecutive period with Temp > ``warm_temp_threshold_c``,
    run ``pond_water_balance`` on that slice only (starting from ``initial_depth_m``). All other
    days use ``initial_depth_m``.
    """
    n = len(meteo_df)
    level_m = np.full(n, float(initial_depth_m), dtype=float)
    dates = pd.to_datetime(meteo_df["Date"])
    years: List[int] = sorted(dates.dt.year.unique().tolist())
    temp = meteo_df["Temp"].to_numpy(dtype=float)
    precip_m = np.asarray(precip_m, dtype=float)
    pet_m = np.asarray(pet_m, dtype=float)

    for year in years:
        ymask = (dates.dt.year == year).to_numpy()
        idx = np.flatnonzero(ymask)
        if idx.size == 0:
            continue
        warm = temp[idx] > float(warm_temp_threshold_c)
        run = longest_consecutive_run_indices(warm)
        if run is None:
            continue
        a, b = run
        g0, g1 = idx[a], idx[b]
        sl = slice(g0, g1 + 1)
        depths = pond_water_balance(
            leaky_params,
            precip_m[sl],
            pet_m[sl],
            initial_depth_m,
        )
        level_m[sl] = depths
    return level_m


def calculate_daily_pet_m(
    temp_c: np.ndarray,
    dates: pd.Series,
    lat_deg: float,
) -> np.ndarray:
    """Hamon PET (m/day) from air temperature (°C) and latitude."""
    day_of_year = pd.to_datetime(dates).dt.dayofyear.to_numpy()
    gamma = 2 * np.pi / 365 * (day_of_year - 1)
    declination = (
        0.006918
        - 0.399912 * np.cos(gamma)
        + 0.070257 * np.sin(gamma)
        - 0.006758 * np.cos(2 * gamma)
        + 0.000907 * np.sin(2 * gamma)
    )
    lat_rad = np.radians(lat_deg)
    cos_hr_angle = -np.tan(lat_rad) * np.tan(declination)
    cos_hr_angle = np.clip(cos_hr_angle, -1.0, 1.0)
    hr_angle = np.arccos(cos_hr_angle)
    daylight_hours = 2 * np.degrees(hr_angle) / 15

    temp_c = np.asarray(temp_c, dtype=float)
    es = 0.6108 * np.exp((17.27 * temp_c) / (temp_c + 237.3))
    pet_mm = (
        2.1 * (daylight_hours / 12.0) ** 2 * (es / (temp_c + 273.3)) * 100.0
    )
    return pet_mm / 1000.0


def load_meteo_dataframe(meteo_path: Union[str, Path]) -> pd.DataFrame:
    """Load meteo CSV; add Date, Temp (°C), Precip (m/day)."""
    df = pd.read_csv(meteo_path)
    df = df.copy()
    df["Date"] = pd.to_datetime(df[["Year", "Month", "Day"]])
    # LAKE meteo: Temp in K -> °C (match leaky_bucket.ipynb: subtract 273)
    df["Temp"] = df["Temp"] - 273.0
    # Precip in file is m/s -> m/day
    df["Precip"] = df["Precip"].astype(float) * 86400.0
    df.loc[df["Precip"] < 0, "Precip"] = 0.0
    return df


@dataclass
class InflowStatic:
    """Static tributary fields (tribheat=2 style); use -999 for unused tracers."""

    width: float = 0.5
    u: float = 0.1
    depth: float = 2.5
    sal: float = -999.0
    ux: float = -999.0
    uy: float = -999.0
    doc: float = -999.0
    poc: float = -999.0
    dic: float = -999.0
    ch4: float = -999.0
    o2: float = -999.0


@dataclass
class OutflowStatic:
    depth0: float = 0.01
    depth: float = 2.5


# Default leaky-bucket parameters from TKL524 calibration in leaky_bucket.ipynb (illustrative).
DEFAULT_LEAKY_PARAMS: Tuple[float, float, float] = (
    5.366202951950946, #3.7780017655431504,
    0.02005233391535097, #0.012612336640402264,
    1.8491948678905278e-07 #2.4629676560539345e-08,
)


def build_inflow_outflow_dataframes(
    meteo_df: pd.DataFrame,
    initial_depth_m: float,
    *,
    latitude_deg: float,
    leaky_params: Tuple[float, float, float] = DEFAULT_LEAKY_PARAMS,
    warm_temp_threshold_c: float = 5.0,
    inflow_static: Optional[InflowStatic] = None,
    outflow_static: Optional[OutflowStatic] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return (inflows_df, outflows_df) ready to write without headers.

    ``meteo_df`` must contain columns Date, Temp (°C), Precip (m/day).

    Leaky-bucket ``level`` uses the longest consecutive warm period per year
    (Temp > ``warm_temp_threshold_c``); other days hold ``initial_depth_m``.
    ``level`` is then smoothed with a centered 5-day rolling mean.
    """
    inflow_static = inflow_static or InflowStatic()
    outflow_static = outflow_static or OutflowStatic()

    meteo_df = meteo_df.sort_values("Date").reset_index(drop=True)
    n = len(meteo_df)

    pet_m = calculate_daily_pet_m(
        meteo_df["Temp"].to_numpy(),
        meteo_df["Date"],
        latitude_deg,
    )
    precip_m = meteo_df["Precip"].to_numpy(dtype=float)
    level_m = level_from_seasonal_leaky_bucket(
        meteo_df,
        precip_m,
        pet_m,
        leaky_params,
        initial_depth_m,
        warm_temp_threshold_c=warm_temp_threshold_c,
    )

    date_str = meteo_df["Date"].dt.strftime("%Y%m%d")

    inflows = pd.DataFrame(
        {
            "Date": date_str,
            "U": inflow_static.u,
            "depth": inflow_static.depth,
            "width": inflow_static.width,
            "temp": meteo_df["Temp"].astype(float),
            "sal": inflow_static.sal,
            "Ux": inflow_static.ux,
            "Uy": inflow_static.uy,
            "DOC": inflow_static.doc,
            "POC": inflow_static.poc,
            "DIC": inflow_static.dic,
            "CH4": inflow_static.ch4,
            "O2": inflow_static.o2,
        }
    )

    outflows = pd.DataFrame(
        {
            "Date": date_str,
            "U": inflow_static.u,
            "depth0": outflow_static.depth0,
            "depth": outflow_static.depth,
            "width": inflow_static.width,
            "level": level_m,
        }
    )
    outflows["level"] = outflows["level"].rolling(
        window=5, min_periods=1, center=True
    ).mean()

    # Match notebook: fill NaN with -999 for inflow tracers (temp is always set)
    inflows = inflows.fillna(-999.0)

    assert len(inflows) == n and len(outflows) == n
    return inflows, outflows


def generate_inflow_outflow(
    meteo_path: Union[str, Path],
    driver_path: Union[str, Path],
    inflow_out_path: Union[str, Path],
    outflow_out_path: Union[str, Path],
    *,
    latitude_deg: Optional[float] = None,
    leaky_params: Tuple[float, float, float] = DEFAULT_LEAKY_PARAMS,
    warm_temp_threshold_c: float = 5.0,
    inflow_static: Optional[InflowStatic] = None,
    outflow_static: Optional[OutflowStatic] = None,
    show_plots: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main API: read meteo + driver, write inflow/outflow files.

    If ``latitude_deg`` is None, ``phi`` is read from the driver file.
    Initial depth is ``h10`` from the driver.

    If ``show_plots`` is True, saves a leaky-bucket timeseries PNG next to this script
    (same directory as ``generate_inflow_outflow.py``) and prints its path.

    ``warm_temp_threshold_c`` (default 5): per-year longest consecutive days with Temp above
    this value get the leaky bucket; all other days use ``h10``.
    """
    meteo_df = load_meteo_dataframe(meteo_path)
    h0 = read_initial_lake_depth_m(driver_path)
    lat = latitude_deg if latitude_deg is not None else read_latitude_deg(driver_path)

    inflows, outflows = build_inflow_outflow_dataframes(
        meteo_df,
        h0,
        latitude_deg=lat,
        leaky_params=leaky_params,
        warm_temp_threshold_c=warm_temp_threshold_c,
        inflow_static=inflow_static,
        outflow_static=outflow_static,
    )

    inflow_out_path = Path(inflow_out_path)
    outflow_out_path = Path(outflow_out_path)
    inflow_out_path.parent.mkdir(parents=True, exist_ok=True)
    outflow_out_path.parent.mkdir(parents=True, exist_ok=True)

    # Same as gen_in_outflows_file.ipynb: no index, no header
    inflows.to_csv(inflow_out_path, index=False, header=False)
    outflows.to_csv(outflow_out_path, index=False, header=False)

    if show_plots:
        mdf = meteo_df.sort_values("Date").reset_index(drop=True)
        saved = save_leaky_bucket_timeseries_plot(
            mdf["Date"],
            mdf["Precip"].to_numpy(),
            mdf["Temp"].to_numpy(),
            outflows["level"].to_numpy(),
            title=str(Path(meteo_path).stem),
        )
        print(f"Saved plot {saved}")

    return inflows, outflows


def _parse_leaky_params(s: str) -> Tuple[float, float, float]:
    parts = [float(x) for x in s.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            "leaky_params must be three comma-separated floats: c_ratio,seepage_m_day,et_scale"
        )
    return (parts[0], parts[1], parts[2])


def main() -> None:
    p = argparse.ArgumentParser(
        description=(
            "Generate LAKE inflow/outflow files from meteo + driver. "
            "Outflow level = h10 outside each year's longest warm spell (see --warm-threshold); "
            "leaky bucket on that spell only, reset to h10 each year."
        )
    )
    p.add_argument("meteo", type=Path, help="Meteo CSV (e.g. Situkuyuk.dat)")
    p.add_argument("driver", type=Path, help="LAKE driver file (e.g. Situkuyuk_driver.dat)")
    p.add_argument("inflow_out", type=Path, help="Output inflows path")
    p.add_argument("outflow_out", type=Path, help="Output outflows path")
    p.add_argument(
        "--lat",
        type=float,
        default=None,
        help="Latitude (deg N) for Hamon PET; default: read ``phi`` from driver",
    )
    p.add_argument(
        "--leaky-params",
        type=_parse_leaky_params,
        default=None,
        help=(
            "Leaky-bucket parameters c_ratio,seepage_m_day,et_scale "
            f"(default: TKL524-calibrated {DEFAULT_LEAKY_PARAMS})"
        ),
    )
    p.add_argument(
        "--warm-threshold",
        type=float,
        default=5.0,
        metavar="C",
        help=(
            "Leaky bucket runs only on each year's longest consecutive period with "
            "air temp (°C) above this value; other days use h10 (default: 5)"
        ),
    )
    p.add_argument("--width", type=float, default=InflowStatic.width)
    p.add_argument("--u", type=float, default=InflowStatic.u, dest="u_vel")
    p.add_argument("--inflow-depth", type=float, default=InflowStatic.depth)
    p.add_argument("--sal", type=float, default=InflowStatic.sal)
    p.add_argument("--ux", type=float, default=InflowStatic.ux)
    p.add_argument("--uy", type=float, default=InflowStatic.uy)
    p.add_argument("--doc", type=float, default=InflowStatic.doc)
    p.add_argument("--poc", type=float, default=InflowStatic.poc)
    p.add_argument("--dic", type=float, default=InflowStatic.dic)
    p.add_argument("--ch4", type=float, default=InflowStatic.ch4)
    p.add_argument("--o2", type=float, default=InflowStatic.o2)
    p.add_argument("--depth0", type=float, default=OutflowStatic.depth0)
    p.add_argument("--outflow-depth", type=float, default=OutflowStatic.depth)
    p.add_argument(
        "--plot",
        action="store_true",
        help=(
            "Save leaky-bucket timeseries PNG next to generate_inflow_outflow.py "
            "(<meteo_stem>_leaky_bucket_timeseries.png); no display required."
        ),
    )

    args = p.parse_args()
    leaky = args.leaky_params if args.leaky_params is not None else DEFAULT_LEAKY_PARAMS

    istatic = InflowStatic(
        width=args.width,
        u=args.u_vel,
        depth=args.inflow_depth,
        sal=args.sal,
        ux=args.ux,
        uy=args.uy,
        doc=args.doc,
        poc=args.poc,
        dic=args.dic,
        ch4=args.ch4,
        o2=args.o2,
    )
    ostatic = OutflowStatic(depth0=args.depth0, depth=args.outflow_depth)

    generate_inflow_outflow(
        args.meteo,
        args.driver,
        args.inflow_out,
        args.outflow_out,
        latitude_deg=args.lat,
        leaky_params=leaky,
        warm_temp_threshold_c=args.warm_threshold,
        inflow_static=istatic,
        outflow_static=ostatic,
        show_plots=args.plot,
    )
    print(f"Wrote {args.inflow_out} and {args.outflow_out}")


if __name__ == "__main__":
    main()

"""
examples

python generate_inflow_outflow.py \
  /home/amullen/LAKE_NSF/meteo/Situkuyuk.dat \
  /home/amullen/LAKE_NSF/setup/Situkuyuk_driver.dat \
  Situkuyuk_inflows.dat Situkuyuk_outflows.dat \
  --plot


python generate_inflow_outflow.py \
  /home/amullen/LAKE_NSF/meteo/Situkuyuk.dat \
  /home/amullen/LAKE_NSF/setup/Situkuyuk_driver.dat \
  Situkuyuk_inflows.dat Situkuyuk_outflows.dat \
  --width 1.0 \
  --u 0.1 \
  --plot

python generate_inflow_outflow.py \
  /home/amullen/LAKE_NSF/meteo/Situkuyuk.dat \
  /home/amullen/LAKE_NSF/setup/Situkuyuk_driver.dat \
  /home/amullen/LAKE_NSF/meteo/Situkuyuk_dynamic_inflows.dat \
  /home/amullen/LAKE_NSF/meteo/Situkuyuk_dynamic_outflows.dat \
  --width 1.0 \
  --u 0.1 \
  --plot
"""
