#!/usr/bin/env python3
"""
Build LAKE-model tributary inflow / outflow ASCII files from a meteo CSV and a driver file.

- Meteo file: same format as LAKE_NSF/meteo/*.dat (Year, Month, Day, Temp in K, Precip in m/s, ...).
- Driver file: LAKE driver; initial liquid water depth is read from keyword ``h10`` (m).
- Inflow temperature follows air temperature (meteo Temp converted to °C).
- Other inflow fields are static scalars (optional arguments).
- Outflow ``level`` (m) is from the leaky-bucket water balance (leaky_bucket.ipynb), using
  Hamon PET from air temperature (same notebook).

Column order matches gen_in_outflows_file.ipynb (no CSV headers):
  Inflow:  Date, U, depth, width, temp, sal, Ux, Uy, DOC, POC, DIC, CH4, O2
  Outflow: Date, U, depth0, depth, width, level
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np
import pandas as pd


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
    """
    c_ratio, seepage_m_day, et_scale = params
    depths = np.zeros(len(precip_m), dtype=float)
    current = float(initial_depth_m)
    for t in range(len(precip_m)):
        inflow = precip_m[t] * (1.0 + c_ratio)
        outflow = seepage_m_day + pet_m[t] * et_scale
        current = max(0.0, current + inflow - outflow)
        depths[t] = current
    return depths


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
    inflow_static: Optional[InflowStatic] = None,
    outflow_static: Optional[OutflowStatic] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return (inflows_df, outflows_df) ready to write without headers.

    ``meteo_df`` must contain columns Date, Temp (°C), Precip (m/day).
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
    level_m = pond_water_balance(leaky_params, precip_m, pet_m, initial_depth_m)

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
    inflow_static: Optional[InflowStatic] = None,
    outflow_static: Optional[OutflowStatic] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main API: read meteo + driver, write inflow/outflow files.

    If ``latitude_deg`` is None, ``phi`` is read from the driver file.
    Initial depth is ``h10`` from the driver.
    """
    meteo_df = load_meteo_dataframe(meteo_path)
    h0 = read_initial_lake_depth_m(driver_path)
    lat = latitude_deg if latitude_deg is not None else read_latitude_deg(driver_path)

    inflows, outflows = build_inflow_outflow_dataframes(
        meteo_df,
        h0,
        latitude_deg=lat,
        leaky_params=leaky_params,
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
        description="Generate LAKE inflow/outflow files from meteo + driver (leaky-bucket level)."
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
        inflow_static=istatic,
        outflow_static=ostatic,
    )
    print(f"Wrote {args.inflow_out} and {args.outflow_out}")


if __name__ == "__main__":
    main()

"""
examples

python generate_inflow_outflow.py \
  /home/amullen/LAKE_NSF/meteo/Situkuyuk.dat \
  /home/amullen/LAKE_NSF/setup/Situkuyuk_driver.dat \
  Situkuyuk_inflows.dat Situkuyuk_outflows.dat


python generate_inflow_outflow.py \
  /home/amullen/LAKE_NSF/meteo/Situkuyuk.dat \
  /home/amullen/LAKE_NSF/setup/Situkuyuk_driver.dat \
  Situkuyuk_inflows.dat Situkuyuk_outflows.dat \
  --width 1.0
  --u 0.1
"""