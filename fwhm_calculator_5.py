from matplotlib.widgets import RectangleSelector
from scipy import interpolate
from scipy import optimize
from scipy.ndimage import affine_transform
from scipy.ndimage import rotate
from skimage.measure import label, regionprops
import matplotlib.pyplot as plt
import numpy as np
import os
import tifffile as tiff

focalSpotImagePath = r"C:\Users\Andy\OneDrive - The Ohio State University\Lab Files\GRAY Laser\Regen\Pump spots\2026-02-22\right 34.tif"
bkgPath = r"C:\Users\Andy\OneDrive - The Ohio State University\Lab Files\2026-01-31 Quantizing Nanolaminates HfO2-SiO2\focal spot\fs_Merge_3.tiff"

bkg = False
subtract_corner = False
roiSelect = False
bkgSelect = True

threshold = 1e-3
px_per_um = 1.0/6.0
lw = 10
power = 4.9  #Watts
rr = 500.0 #Hz


roi = {"x0": None, "x1": None, "y0": None, "y1": None}

def peak_fluence(power, rr, fwhm_x, fwhm_y):
    return (4.0 * (power / rr) * np.log(2.0)) / (np.pi * fwhm_x * fwhm_y)

def px_to_um(px):
    return px*(1.0/px_per_um)

def um_to_cm(um):
    return um * 1e-4

def normalize_to_1(data):
    return data / np.max(data)

def corner_subtract(data, corner_size=100):
    corner0 = data[:corner_size,:corner_size]
    corner1 = data[:corner_size,-corner_size:]
    corner2 = data[-corner_size:,:corner_size]
    corner3 = data[-corner_size:,-corner_size:]
    return data - (np.mean(corner0 + corner1 + corner2 + corner3) / 4.0)

def select_roi(data, title="Drag ROI, then close window", minspanx=5, minspany=5, button=1):
    roi = {"x0": None, "x1": None, "y0": None, "y1": None}

    def onselect(eclick, erelease):
        x0, y0 = int(eclick.xdata), int(eclick.ydata)
        x1, y1 = int(erelease.xdata), int(erelease.ydata)
        roi["x0"], roi["x1"] = sorted([x0, x1])
        roi["y0"], roi["y1"] = sorted([y0, y1])

    fig, ax = plt.subplots()
    ax.imshow(data, cmap='inferno')
    rs = RectangleSelector(ax, onselect, useblit=True, interactive=True,
                        button=[1], minspanx=5, minspany=5, spancoords="pixels")
    plt.title(title)
    plt.show()

    return roi["x0"], roi["x1"], roi["y0"], roi["y1"]
    
def lineouts(data, lw=0, cy=None, cx=None):
    #lw is the amount of pixels on either side of the main lineout
    #TODO: add in option for averaging lineouts > 1 pixel wide
    if ((cy is None) and (cx is None)):
        y0, x0 = np.unravel_index(np.argmax(data), data.shape)
    else:
        y0=int(np.round(cy))
        x0=int(np.round(cx))

    if (lw > 0):
        x_lineout = np.mean(data[:,x0-lw:x0+lw], axis=1)
        y_lineout = np.mean(data[y0-lw:y0+lw,:], axis=0)
    else:
        x_lineout = data[:,x0]
        y_lineout = data[y0,:]
    return x_lineout, y_lineout

def calc_fwhm(data):
    pixels = np.arange(len(data))

    cubic_spline = interpolate.CubicSpline(pixels, data)
    
    res = optimize.minimize_scalar(lambda x: -cubic_spline(x), bounds=(np.min(pixels), np.max(pixels)))
    half_max = cubic_spline(res.x)/2.0
    half_max = np.max(data)/2.0
    gray_arg_max = np.argmax(data)
    
    
    pixels1 = pixels[:gray_arg_max]
    pixels2 = pixels[gray_arg_max:]
    
    grays1 = data[:gray_arg_max]
    grays2 = data[gray_arg_max:]
    
    f1 = interpolate.interp1d(grays1, pixels1)
    f2 = interpolate.interp1d(grays2, pixels2)
    
    #return(abs(f1(half_max)-f2(half_max)))
    return f1(half_max), f2(half_max)

def rotate_about_center(data, theta, cx, cy, order=3):
    """
    Rotate image I by angle theta (radians) about (cx, cy) in pixel coords.
    Returns rotated image with same shape.
    """
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s],
                  [s,  c]])  # forward rotation in (x,y)

    # scipy.ndimage.affine_transform maps output coords -> input coords
    # so use inverse rotation:
    A = R.T

    center = np.array([cx, cy])  # (x, y)
    offset = center - A @ center

    # affine_transform expects matrix in (row,col) = (y,x) order
    A_rc = np.array([[A[1,1], A[1,0]],
                     [A[0,1], A[0,0]]])
    offset_rc = np.array([offset[1], offset[0]])

    return affine_transform(data, A_rc, offset=offset_rc, order=order, mode="constant", cval=0.0)

if bkg:
    data = tiff.imread(focalSpotImagePath).astype(float) - tiff.imread(bkgPath).astype(float)
else:
    data = tiff.imread(focalSpotImagePath).astype(float)
    if subtract_corner:
        data = corner_subtract(data)
if bkgSelect:
    bkg_x0, bkg_x1, bkg_y0, bkg_y1 = select_roi(data, title='Drag background, then close window')
    data -= np.mean(data[bkg_y0:bkg_y1, bkg_x0:bkg_x1])
    data = np.clip(data, a_min=0, a_max=None)

if roiSelect:
    signal_x0, signal_x1, signal_y0, signal_y1 = select_roi(data, title='Drag signal, then close window')
    data = data[signal_y0:signal_y1, signal_x0:signal_x1]

mask0 = data > (threshold * np.max(data))
lbl = label(mask0)

py, px = np.unravel_index(np.argmax(data), data.shape)
lab_peak = lbl[py, px]
mask = (lbl == lab_peak)

props = regionprops(mask.astype(int), intensity_image=data)[0]

cy, cx = props.centroid_weighted
theta = props.orientation
theta_deg = np.degrees(theta)

mu20 = props.moments_weighted_central[2, 0]
mu02 = props.moments_weighted_central[0, 2]
mu11 = props.moments_weighted_central[1, 1]

M00 = props.weighted_moments[0, 0]
Sigma = np.array([[mu20, mu11], [mu11, mu02]]) / M00

evals, evecs = np.linalg.eigh(Sigma)
w = 2.0 * np.sqrt(evals)
fwhm_moments_px = w * np.sqrt(2.0 * np.log(2.0))
fwhm_moments_um = px_to_um(fwhm_moments_px)
# Coordinates of the max intensity pixel
# max_idx = np.argmax(data)
# y0_init, x0_init = np.unravel_index(max_idx, data.shape)
# ny_full, nx_full = data.shape

# y_min = max(y0_init - crop_hw, 0)
# y_max = min(y0_init + crop_hw, ny_full)
# x_min = max(x0_init - crop_hw, 0)
# x_max = min(x0_init + crop_hw, nx_full)

# plt.imshow(data)
# plt.axhline(y=y0_init, c='C1')
# plt.axvline(x=x0_init, c='C1')
# plt.show()

#data = data[y_min:y_max, x_min:x_max]

#data_rot = rotate(data, angle=theta_deg, reshape=False)
data_rot = rotate_about_center(data, theta, cx, cy)
# plt.imshow(data_rot, cmap='inferno')
# plt.title('rotated')
# plt.show()

x_lineout, y_lineout = lineouts(data_rot, lw=lw, cx=cx, cy=cy)

x_lineout = normalize_to_1(x_lineout)
y_lineout = normalize_to_1(y_lineout)

# plt.imshow(data)
# plt.title('unrotated')
# plt.show()



x1, x0 = calc_fwhm(x_lineout)
y1, y0 = calc_fwhm(y_lineout)

fwhm_x_px = abs(x1-x0)
fwhm_y_px = abs(y1-y0)

fwhm_x_um = px_to_um(fwhm_x_px)
fwhm_y_um = px_to_um(fwhm_y_px)

fwhm_multiplier = 1.5
plot_lower_x = np.argmax(x_lineout) - fwhm_multiplier * fwhm_x_px
plot_upper_x = np.argmax(x_lineout) + fwhm_multiplier * fwhm_x_px
plot_lower_y = np.argmax(y_lineout) - fwhm_multiplier * fwhm_y_px
plot_upper_y = np.argmax(y_lineout) + fwhm_multiplier * fwhm_y_px

plt.figure(figsize=(12,4))
plt.subplot(1, 3, 1)
plt.imshow(data_rot, cmap='inferno')
plt.xlim((plot_lower_y, plot_upper_y))
plt.ylim((plot_lower_x, plot_upper_x))

plt.subplot(1, 3, 2)
plt.plot(np.arange(len(x_lineout)),x_lineout)
plt.xlabel('pixel index')
plt.ylabel('intensity')
plt.axhline(y=np.max(x_lineout)/2.0)
plt.axvline(x=x0)
plt.axvline(x=x1)
plt.xlim((plot_lower_x, plot_upper_x))
plt.title(f'x lineout averaged over {2*lw+1} pixels')

plt.subplot(1, 3, 3)

plt.plot(np.arange(len(y_lineout)),y_lineout)
plt.xlabel('pixel index')
plt.ylabel('intensity')
plt.axhline(y=np.max(y_lineout)/2.0)
plt.axvline(x=y0)
plt.axvline(x=y1)
plt.xlim((plot_lower_y, plot_upper_y))
plt.title(f'y lineout averaged over {2*lw+1} pixels')
plt.show()

# print("From moments")
# print(f"x FWHM: {fwhm_moments_px[0]} pixels")
# print(f"y FWHM: {fwhm_moments_px[1]} pixels")
# print()
# print(f"x FWHM: {fwhm_moments_um[0]} microns")
# print(f"y FWHM: {fwhm_moments_um[1]} microns")
# print()
print("From Lineouts")
print(f"x FWHM: {fwhm_x_px} pixels")
print(f"y FWHM: {fwhm_y_px} pixels")
print()
print(f"x FWHM: {fwhm_x_um} microns")
print(f"y FWHM: {fwhm_y_um} microns")
print()
print(f"ratio: {np.max(np.array([fwhm_x_px, fwhm_y_px]))/np.min(np.array([fwhm_x_px, fwhm_y_px]))}")
print()
print(f"fluence: {peak_fluence(power, rr, um_to_cm(fwhm_x_um), um_to_cm(fwhm_y_um)):.2f} J/cm^2")