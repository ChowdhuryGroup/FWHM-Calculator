# FWHM-Calculator

A small Python utility to measure the full width at half maximum (FWHM) of a focal spot (or bright region) in a TIFF image. The script locates the brightest object, computes axis-aligned and rotated lineouts about the weighted centroid, fits/interpolates the lineouts to find FWHM values (in pixels and microns), and computes peak fluence for pulsed sources.

This repository contains a single-script workflow useful for analyzing laser focal spots and similar intensity distributions.

---

## Features

- Load TIFF images (single-frame) via tifffile
- Optional background subtraction (image or ROI-based)
- Optional corner background subtraction
- Automatic region detection (threshold + connected component labeling)
- Weighted centroid and orientation via skimage regionprops
- Rotation about the centroid to align the major axes
- Lineouts (with averaging width) and cubic-spline-based FWHM calculation
- FWHM reported in pixels and microns (user-supplied scale)
- Peak fluence calculation given pulse energy/pulse repetition rate
- Visual output: rotated image and lineout plots showing half-maximum and FWHM markers

---

## Quick Start / Requirements

Install the Python dependencies:

```bash
pip install numpy scipy matplotlib scikit-image tifffile
```

(Developed/tested with numpy, scipy, matplotlib, scikit-image and tifffile.)

---

## How to use

1. Put your image file path in the `focalSpotImagePath` variable near the top of the script.
2. (Optional) Put a background image path in `bkgPath`, or enable `bkgSelect` to select a background ROI interactively.
3. Adjust parameters such as `threshold`, `px_per_um`, `lw`, `power`, and `rr`.
4. Run the script from a Python interpreter:

```bash
python fwhm_calculator.py
```

During execution:
- If `bkgSelect = True` the script will prompt you to drag a rectangle in a displayed image to choose a background ROI (click-drag and close the window).
- If `roiSelect = True` the script will prompt you to select the signal ROI similarly.

The script then shows a 3-panel figure:
- Rotated image (zoomed to a few times the FWHM)
- x lineout with half-maximum and FWHM markers
- y lineout with half-maximum and FWHM markers

It also prints FWHM and peak fluence values to the console.

---

## Key parameters (top of script)

- focalSpotImagePath (str): Path to the TIFF focal spot image.
- bkgPath (str): Path to a background TIFF to subtract (if `bkg=True`).
- bkg (bool): If True, subtract `bkgPath` image from `focalSpotImagePath`.
- subtract_corner (bool): If True, subtract mean of image corners as background.
- roiSelect (bool): If True, interactively select a signal ROI to crop before analysis.
- bkgSelect (bool): If True, interactively select a background ROI and subtract its mean.
- threshold (float): Fraction of the maximum (0-1) used to create a mask for region detection (e.g., 1e-3).
- px_per_um (float): Pixel density, pixels per micron (default shown is 1/6.0 -> 6 px/µm).
- lw (int): Half-width (in pixels) for averaging lineouts (averages 2*lw+1 pixels).
- power (float): Average laser power in Watts (used to compute energy per pulse).
- rr (float): Repetition rate in Hz (used to compute energy per pulse).

Units:
- FWHM printed in pixels and microns (uses `px_per_um`).
- Peak fluence printed in J/cm^2 (converts microns to cm inside the calculation).

---

## Important functions (what they do)

- peak_fluence(power, rr, fwhm_x, fwhm_y):
  Computes peak fluence using a Gaussian focal spot model:
  fluence = 4 * (power/rr) * ln(2) / (π * fwhm_x * fwhm_y)

- px_to_um(px) / um_to_cm(um):
  Unit conversion helpers.

- corner_subtract(data, corner_size=100):
  Estimate background from four image corners and subtract it.

- select_roi(data, title=..., minspanx=5, ...):
  Interactive rectangle selection using matplotlib.widgets.RectangleSelector. Returns x0,x1,y0,y1.

- lineouts(data, lw=0, cy=None, cx=None):
  Returns 1D lineouts through the center (or provided centroid) optionally averaged over a small width.

- calc_fwhm(data):
  Given a 1D intensity vector, finds interpolated left/right positions at half-maximum and returns them (uses cubic interpolation and scipy).

- rotate_about_center(data, theta, cx, cy):
  Rotate the image around a given centroid using an affine transform while preserving image shape.

---

## Typical workflow (what the script does)

1. Read image (optionally subtract background image or corner mean).
2. Optional background ROI selection and subtraction.
3. Optional signal ROI selection and crop.
4. Mask image with `data > threshold * max(data)` and label connected components.
5. Find the labeled object containing the global maximum; compute weighted centroid and orientation.
6. Compute central moments to estimate elliptical widths (moments-based FWHM).
7. Rotate image about centroid so the spot major/minor axes align with image axes.
8. Extract averaged lineouts, normalize, interpolate, and compute FWHM by locating half-maximum crossings.
9. Plot results and print numeric outputs including peak fluence.

---

## Output example (console)

```
From Lineouts
x FWHM: 23.4 pixels
y FWHM: 17.8 pixels

x FWHM: 3.9 microns
y FWHM: 3.0 microns

ratio: 1.31

fluence: 0.12 J/cm^2
```

(Values are illustrative — actual values vary with image & parameters.)

---

## Tips & troubleshooting

- If the script finds the wrong object, adjust `threshold` or crop the image with `roiSelect=True`.
- Use `lw` > 0 to average lineouts and reduce noise before FWHM interpolation.
- If the spot extends to the image edge, consider padding/cropping to avoid truncation errors.
- For very noisy images, smooth the image or the lineouts (e.g., with a small Gaussian) prior to FWHM calculation, but be aware smoothing changes apparent FWHM.

---

## Next improvements you might consider

- Add a small CLI (argparse) so users can pass file paths and parameters without editing the script.
- Add unit tests for the FWHM routine using synthetic Gaussian spots.
- Add support for multi-page TIFF stacks or batch processing of folders.
- Save result images (PNG) and a CSV summary of measured values.

---

## Contributing

If you'd like contributions:
- Open an issue describing the feature or bug.
- I can help convert the script into a CLI module or open a PR — tell me which you prefer.

---

## License & authorship

Add a license file if you want to make this open source (MIT/BSD recommended for small scripts). If you want, I can add a simple MIT license file for the repository.

---

Thank you — this README documents the included FWHM measurement script, explains usage and parameters, and suggests practical next steps.
