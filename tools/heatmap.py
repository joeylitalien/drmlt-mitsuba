from __future__ import print_function

"""
A simple script to compute heatmap of accepted stages.
"""

import argparse
import pyexr
import numpy as np
import matplotlib.pyplot as plt

COLOR_MAP = 'plasma'

def stage_map(img, clip, eps=1e-2):
    """Compute ratio of stage acceptances."""

    first_stage = img[:,:,0]
    second_stage = img[:,:,1]
    normal = first_stage + second_stage + eps
    for i in range(3):
        img[:,:,i] /= normal

    return img[:,:,1]

def plot(img, clip, fname):
    """Plot false color heatmap with colorbar legend."""

    plt.rcParams.update({'font.size': 12})
    plt.rcParams.update({'font.family': 'linux biolinum'})

    dpi_scale = 60.0
    figsize = img.shape[0] / dpi_scale, img.shape[1] / dpi_scale

    fig, ax1 = plt.subplots(1, 1, figsize=figsize)
    ax1.set_axis_off()
    ax1.xaxis.set_major_locator(plt.NullLocator())
    ax1.yaxis.set_major_locator(plt.NullLocator())

    img1 = plt.imshow(img, cmap=COLOR_MAP, interpolation='nearest', vmin=clip[0], vmax=clip[1])

    
    plt.tight_layout(h_pad=1)
    fig.savefig(fname, bbox_inches='tight', pad_inches=0.1)


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='Compute delayed rejection acceptance map.')
    parser.add_argument('-t',   '--test', help='test image filename', type=str, required=True)
    parser.add_argument('-eps', '--epsilon', help='epsilon value', type=float, default=1e-2)
    parser.add_argument('-c',   '--clip', help='clipping values for min/max', nargs=2, type=float, default=[0, 1])
    parser.add_argument('-cb',  '--colorbar', help='display colorbar on false error heatmap', action='store_true')

    args = parser.parse_args()

    # Check image format
    if not  args.test.lower().endswith('.exr'):
        raise ValueError('Image must be in OpenEXR format.')

    # Load images and convert to NumPy array
    try:
        test_fp = pyexr.open(args.test)
    except FileNotFoundError:
        print('Could not open file')
    test = np.array(test_fp.get())

    stage_img = stage_map(test, args.clip, args.epsilon)
    fname = 'acceptance-map.png'
    plot(stage_img, args.clip, fname)