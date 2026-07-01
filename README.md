# MANB

MATLAB code for **Ultra-Narrow Bands and High Sensitivity in Moiré Acoustic Waveguides**.

## Repository layout

- `functions/` — shared utilities (transfer-matrix method, band structure, IPR, narrow-band detection, optimization)
- `Moire_TL.m`, `Moire_Band_Narrow.m`, `Moire_Ipr.m`, etc. — main simulation and analysis scripts
- `moire_*_scan_parfor.m` — parameter scans (parallel)

## Requirements

- MATLAB (R2018b or later recommended)
- Parallel Computing Toolbox for `parfor` scan scripts

## Usage

1. Open MATLAB and set the current folder to this repository root.
2. Add `functions` to the path, or run scripts that call `addpath` as needed.
3. Run the desired script (e.g. `Moire_TL.m`).

## Citation

If you use this code, please cite the associated paper on ultra-narrow bands and sensitivity in moiré acoustic waveguides.
