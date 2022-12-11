"""Identify 3D centromere-bounding regions of interest and quantify data within them.

Imports:

    os
    sys
    datetime.date
    numpy
    scipy.ndimage
    skimage.io.imread
    skimage.feature.peak_local_max
    roifile

Classes:

    CentromereQuantifier(numpy arr, numpy arr, string, int):
        Variables:
            condition
            cell_number
            mask_input
            data_input
            dapi_input
            centromere_area_array
            filtered_array
            maxima_list
            roi_array
            bg_array
            final_data
            centromere_radius
            xscale
            yscale
            zscale
            min_centromere_dist
            max_num_maxima
            maxima_std_threshold
            roi_xradius
            roi_yradius
            roi_zradius
            dist_from_edge_xy
            dist_from_edge_z
            bg_xradius
            bg_yradius
            bg_zradius
            bg_std_threshold
        Functions:
            identify_centromere_region() -> numpy array
            identify_centromeres() -> numpy array
            find_maxima() -> list
            draw_labeled_rois() -> numpy array
            draw_labeled_bg() -> numpy array
            quantify() -> numpy array
            save_data() -> None
            save_rois() -> None

Static Functions:

    get_input(string) -> list of couplets of numpy arrays
    rasterize(tuple, list of lists, list) -> numpy array
"""

import os
import sys
from datetime import datetime

import numpy as np
from scipy import ndimage as ndi
from skimage.io import imread
from skimage.feature import peak_local_max
from skimage.filters import threshold_isodata
import roifile


class CentromereQuantifier:
    def __init__(self, condition, cell_number, mask_input, data_input, dapi_input=None):
        self.condition = condition
        self.cell_number = cell_number
        self.mask_input = mask_input
        self.data_input = data_input
        self.dapi_input = dapi_input

        self.centromere_area_array = None
        self.filtered_array = None
        self.maxima_list = None
        self.roi_array = None
        self.bg_array = None
        self.final_data = None

        self.centromere_radius = 150  # nm; half width at half maximum; set empirically
        self.xscale = 147  # nm / pixel
        self.yscale = 147  # nm / pixel
        self.zscale = 500  # nm / slice
        self.min_centromere_dist = 2  # pixels; Chebyshev distance
        self.max_num_maxima = 38  # starting from the most intense local maximum
        self.maxima_std_threshold = 0  # minimum number of standard deviations above the mean to accept a maximum
        self.roi_xradius = 4  # pixels
        self.roi_yradius = 4  # pixels
        self.roi_zradius = 3  # pixels
        self.dist_from_edge_xy = 20  # pixels; minimum distance from edge of image for an ROI
        self.dist_from_edge_z = 2  # slices; minimum distance from edge of image for an ROI
        self.bg_xradius = self.roi_xradius + 1
        self.bg_yradius = self.roi_yradius + 1
        self.bg_zradius = self.roi_zradius + 1
        self.bg_std_threshold = None  # background standard deviations below the mean to exclude a centromere.

    def identify_centromere_region(self):
        """Create a binary mask separating the centromere area from the surrounding background using DAPI image.

        If no DAPI image is provided, this mask is an array of the same shape as the mask array with all values 1.
        """
        if self.dapi_input is None:
            centromere_mask = np.full_like(self.mask_input, 1)
        else:
            dapi_dilated = ndi.grey_dilation(self.dapi_input, size=(3, 16, 16), mode='nearest')
            centromere_mask = dapi_dilated > threshold_isodata(dapi_dilated)
        self.centromere_area_array = centromere_mask
        return centromere_mask

    def identify_centromeres(self):
        """Use a Difference of Gaussians algorithm to isolate centromeres from the assigned mask array."""
        imagearr = self.mask_input
        xbase = self.centromere_radius / self.xscale
        ybase = self.centromere_radius / self.yscale
        zbase = self.centromere_radius / self.zscale
        convolutions = []
        for exp in [0, 1, 2]:
            xsigma = xbase * np.sqrt(2) ** exp
            ysigma = ybase * np.sqrt(2) ** exp
            zsigma = zbase * np.sqrt(2) ** exp
            convolutions.append(ndi.gaussian_filter(imagearr,
                                                    [3 * zsigma, xsigma, ysigma],  # SHOULD X AND Y BE SWITCHED?
                                                    output=float,
                                                    mode='nearest'))
        product = (imagearr-convolutions[0]) * (convolutions[0]-convolutions[1]) * (convolutions[1]-convolutions[2])
        self.filtered_array = product
        return product

    def find_maxima(self):
        """Create a list of indices of the local maxima of a given image array."""
        imagearr = self.filtered_array * self.centromere_area_array
        threshold = np.mean(imagearr)+self.maxima_std_threshold*np.std(imagearr)
        maxima = peak_local_max(imagearr,
                                min_distance=self.min_centromere_dist,
                                num_peaks=self.max_num_maxima,
                                threshold_abs=threshold,
                                exclude_border=(self.dist_from_edge_z, self.dist_from_edge_xy, self.dist_from_edge_xy))
        self.maxima_list = maxima
        return maxima

    def draw_labeled_rois(self):
        shape = self.filtered_array.shape
        centers_of_mass = self.maxima_list
        roi_radii = np.array([self.roi_zradius, self.roi_xradius, self.roi_yradius])
        rois = rasterize(shape, centers_of_mass, roi_radii)
        self.roi_array = rois
        return rois

    def draw_labeled_bg(self):
        """Use rasterize to define background region, removing pixels where the background overlaps with ROIs."""
        shape = self.filtered_array.shape
        centers_of_mass = self.maxima_list
        bg_radii = np.array([self.bg_zradius, self.bg_xradius, self.bg_yradius])
        bgs = rasterize(shape, centers_of_mass, bg_radii)
        bgs[np.where(self.roi_array != 0)] = 0
        self.bg_array = bgs
        return bgs

    def quantify(self):
        def quantify_region(vals):
            average = np.mean(vals)
            npixels = len(vals)
            return [average, npixels]

        mask = self.roi_array
        bg_mask = self.bg_array
        data = self.data_input

        num_centromeres = np.amax(mask)
        label_vals = np.arange(1, num_centromeres + 1)
        # Could remove some of these transpositions, probably.
        info = np.transpose(np.broadcast_to([self.condition, self.cell_number], (num_centromeres, 2)))
        centromere_data = np.transpose(np.asarray([quantify_region(data[mask == value]) for value in label_vals]))
        bg_data = np.transpose(np.asarray([quantify_region(data[bg_mask == value]) for value in label_vals]))
        difference = np.array(centromere_data[0] - bg_data[0])[np.newaxis]
        normalized = np.array(difference[0] / bg_data[0])[np.newaxis]
        raw = np.concatenate((info, centromere_data, bg_data, difference, normalized), axis=0)
        # Take out outliers with extremely low background values, if a threshold is given:
        if self.bg_std_threshold is None:
            final = np.transpose(raw)
        else:
            bgs = np.asarray(raw[4], dtype=float)
            final = [x for x in np.transpose(raw) if float(x[4]) >= np.mean(bgs)-self.bg_std_threshold*np.std(bgs)]
        self.final_data = final
        return final

    def run(self):
        self.identify_centromere_region()
        self.identify_centromeres()
        self.find_maxima()
        self.draw_labeled_rois()
        self.draw_labeled_bg()
        self.quantify()

    def save_data(self, output_file):
        with open(output_file, 'ta') as file:
            np.savetxt(file, self.final_data, fmt='%s,%s,%.7s,%s,%.7s,%s,%.7s,%.5s', delimiter=',')

    def save_rois(self, output_dir):
        """Create an ImageJ ROI zip file from the maxima."""
        condition = self.condition
        cell_number = self.cell_number
        maxima_list = self.maxima_list

        def make_roi_points(coordinates):
            rois = []
            for ind, coord in enumerate(coordinates):
                current = roifile.ImagejRoi()
                current.version = 228
                # Add one because ImageJ slice indexing starts at 1:
                current.position = coord[0] + 1
                # Switch x and y because that's the way it is:
                current.left = coord[2]
                current.top = coord[1]

                current.right = current.left + 1
                current.bottom = current.top + 1
                current.n_coordinates = 1
                current.stroke_width = 2
                current.roitype = roifile.ROI_TYPE.POINT
                current.name = 'ROI ' + str(ind + 1) + ':[' + str(current.position) + ',' + str(
                    current.left) + ',' + str(current.top) + ']'
                rois.append(current)
            return rois

        roi_zip = make_roi_points(maxima_list)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        output_path = os.path.join(output_dir, condition + '_s' + str(cell_number) + '_ROI.zip')
        roifile.roiwrite(output_path, roi_zip, mode='w')


def get_input(directory):
    """ Load a tif file or files into a couplet of lists of numpy arrays, of (mask_images, data_images).

    :param directory: The directory to retrieve from.
    :return: A list of numpy array corresponding to each tif file.
    """
    full_path = os.path.join(os.getcwd(), directory)
    if not os.path.isdir(full_path):
        print('Invalid directory path. Exiting.')
        sys.exit(1)
    mask_regex = input('Mask file pattern with an * in place of cell number: ')
    data_regex = input('Data file pattern with an * in place of cell number: ')
    dapi_regex = input('Dapi file pattern with an * in place of cell number (optional): ')
    if '*' not in mask_regex or '*' not in data_regex or (dapi_regex != '' and '*' not in dapi_regex):
        print('File patterns must contain an *. Exiting.')
        sys.exit(2)

    cellnum = 1
    mask_list = []
    data_list = []
    dapi_list = []
    while os.path.isfile(os.path.join(full_path, mask_regex.replace('*', str(cellnum)))) and\
          os.path.isfile(os.path.join(full_path, data_regex.replace('*', str(cellnum)))):
        mask_file = os.path.join(full_path, mask_regex.replace('*', str(cellnum)))
        mask_list.append(imread(mask_file))
        data_file = os.path.join(full_path, data_regex.replace('*', str(cellnum)))
        data_list.append(imread(data_file))
        if dapi_regex == '':
            dapi_list.append(None)
        else:
            dapi_file = os.path.join(full_path, dapi_regex.replace('*', str(cellnum)))
            dapi_list.append(imread(dapi_file))
        cellnum += 1
    print(f'{cellnum - 1} files found.')
    return mask_list, data_list, dapi_list


def rasterize(shape, coordinate_list, radii):
    """Converts a list of coordinates into a numpy array of pixelated ellipsoids of given radii.

    From https://stackoverflow.com/a/69444906/17045291.
    :param shape: A tuple of array dimensions for the output array
    :param coordinate_list: a list of N-item lists with the coordinates in order of [z, x, y].
    :param radii: an N-item list storing the radii as [z_radius, x_radius, y_radius].
    :return: An N-dimensional numpy array of labeled ellipsoids.
    """
    sh = shape
    out = np.zeros(sh, int)
    aux = np.zeros(sh)
    radii = radii
    for j, com in enumerate(coordinate_list, start=1):
        bboxl = np.floor(com - radii).clip(0, None).astype(int)
        bboxh = (np.ceil(com + radii) + 1).clip(None, sh).astype(int)
        roi = out[tuple(map(slice, bboxl, bboxh))]
        roiaux = aux[tuple(map(slice, bboxl, bboxh))]
        logrid = *map(np.square, np.ogrid[tuple(
            map(slice, (bboxl - com) / radii, (bboxh - com - 1) / radii, 1j * (bboxh - bboxl)))]),
        dst = (1 - sum(logrid)).clip(0, None)
        mask = dst > roiaux
        roi[mask] = j
        np.copyto(roiaux, dst, where=mask)
    return out


if __name__ == "__main__":
    data_dir = input("Path to directory containing mask and data images: ")
    input_triple = get_input(data_dir)
    mask_image_list, data_image_list, dapi_image_list = input_triple

    cond = input('Condition Name: ')
    now = datetime.now().strftime('%Y-%m-%d_%H-%M')
    output_filename = now + '_' + cond
    output_csv = output_filename if output_filename.endswith('.csv') else output_filename + '.csv'
    roi_dir = os.path.join(now + '_' + cond + '_ROIs')
    np.savetxt(output_csv, [],
               header='Condition,Cell,'
                      'Centromere Mean,Centromere Pixels,'
                      'Background Mean,Background Pixels,'
                      'Centromere-Background,Normalized', comments='')

    for i, images in enumerate(zip(mask_image_list, data_image_list, dapi_image_list), start=1):
        mask_image, data_image, dapi_image = images
        print('Quantifying cell ' + str(i) + '...')
        quant = CentromereQuantifier(cond, i, mask_image, data_image, dapi_image)
        quant.run()
        quant.save_data(output_csv)
        quant.save_rois(roi_dir)
    print('Quantification complete.')
    sys.exit(0)


# def centrocalc(outfile, condition_name, cell_number, mask_image, data_image):
#     """ Saves a .csv file of the data from a data image, using a mask image to identify centromeres.
#
#     :param output: A string with the output file name.
#     :param cell_number: An integer with the cell number.
#     :param mask_image: A numpy array loaded from the mask image.
#     :param data_image: A numpy array loaded from the data image.
#     :return: None.
#     """
#     maxima = find_peaks(gaussian_filter(mask_image))
#     labeled_mask = draw_labeled_rois(mask_image, maxima)
#     data = quantify(labeled_mask, data_image)
#     # Add condition name and cell number to each row:
#     info = [condition_name, cell_number]
#     appended = [info + row for row in data]
#     # Take out outliers with background values more than 3 standard deviations less than the mean:
#     bgs = np.asarray(np.transpose(appended)[3], dtype=int)
#     bg_average = np.mean(bgs)
#     bg_std = np.std(bgs)
#     trimmed = [row for row in appended if int(row[3]) >= bg_average - 3 * bg_std]
#     # Save data and regions of interest:
#     save_to(outfile, trimmed)
#     save_rois(condition_name, cell_number, maxima)
