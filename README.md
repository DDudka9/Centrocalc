# 3D Centromere Quantification

Author: R. Brian Akins

Last updated: 8 Nov 2021.


## Set up
### Python Environment
Development was performed using Python 3.9.
Required third-party packages are listed in requirements.txt. These are: 
numpy, scipy, scikit-image, and roifile. All of these can be downloaded from Anaconda.

To download the lastest version of Miniconda, visit: https://docs.conda.io/en/latest/miniconda.html.

To create an environment with these modules, use the terminal to run:
```
conda env create --name centrocalc --file centrocalc.yml
```
Alternatively:
```
conda create --name centrocalc numpy scipy scikit-image roifile python=3.9
```
Once an environment has been created, it can be activated with:
```
conda activate centrocalc
```
### File preparation
This program automatically iterates through image files to quantify multiple cells.
It requires a standardized file naming system where the only difference between files of the same condition
is the cell number. It is case-sensitive. For example:
```
CENPO_GFP_s1.TIF
CENPO_GFP_s2.TIF
CENPO_GFP_s3.TIF
...
```
will work. However,
```
CENPO_GFP_s1.TIF
CENPO_GFP_s2.tif     # Different case in file extension
Cenpo_GFP_s3.TIF     # Different case in file name
CENPO_GFP_cell4.TIF  # Different file name
```
would all be considered different conditions and will not work.

These files should all be stored in the same folder/directory.
## Running quantifications
After activating the conda environment, navigate to the directory containing the program with ```cd``` 
and run it with:
```
python centrocalc.py
```
If the program runs correctly, it will show a prompt saying 
```
Path to directory containing mask and data images: 
```
Input a relative or absolute path to the directory containing both the images to be used as a mask 
and those to be quantified. Then, at the next prompt,
```
Mask file pattern with an * in place of cell number: 
```
enter a string matching the file name of the images to be used as the mask. From our previous example, we would enter:
```
Mask file pattern with an * in place of cell number: CENPO_GFP_s*.TIF
```
Then, do the same for the data images, which may be the same files as the mask images if desired.

Finally, input a condition name to be associated with these data.

The final output will be:
1. a .csv file with the data.
2. a folder of .zip files containing ImageJ ROIs. These contain the coordinates used for each centromere and 
may be viewed in ImageJ.

## Operation Overview
The program operates in six steps.
### First, segment centromere area
If a DAPI stain image is given, it is used to isolate the centromeres from the background. The dapi image is processed
with a greyscale dilation at a width of 16 pixels to ensure the centromeres are included. Then, a threshold is 
calculated using the ISODATA method such that "returned thresholds are intensities that separate the image into two 
groups of pixels, where the threshold intensity is midway between the mean intensities of these groups." The lowest of 
these thresholds is used to separate the cell area from the image background.
If no dapi image is given, the entire cell is considered.
### Second, isolate centromeres
This step is inspired by Vermolen et al., 2008. A difference of Gaussians algorithm is used to isolate spots of a 
known size from the image. Three 3D gaussian convolutions are performed with widths: CENTROMERE_RADIUS, 
CENTROMERE_RADIUS\*sqrt(2), and CENTROMERE_RADIUS\*2. The output is a product of the differences between these 
gaussians, isolating spots of centromere size.

NB: The point spread function of confocal microscopes is 3 times larger in the Z direction than in the X or Y direction,
and the 3D gaussian is scaled to match.
### Third, find local maxima
A local maxima algorithm is used to identify centromeres. Up to MAX_NUM_MAXIMA spots are chosen, separated by a minimum
of MIN_CENTROMERE_DIST pixels by Chebyshev distance (both orthogonal and diagonal directions are of length 1). 
They must have an absolute value of at least MAXIMA_STD_THRESHOLD standard deviations above the mean
and be at least DIST_FROM_EDGE away from the edges of the image (in the x, y, and z directions independently).
### Fourth, draw regions of interest
Using these local maxima points, 3D ellipsoid ROIs are drawn of a given radius. Background ROIs are also drawn as
larger ellipsoids with the centromere ROI excluded. The function to draw ROIs, rasterize(shape, coordinate_list, radii),
is from a Stack Overflow answer at https://stackoverflow.com/a/69444906/17045291.

Overlapping ROIs are resolved by distance, where each pixel is assigned to the closest maxima. These overlap resolution
distances are non-Euclidean; they are scaled to the radius of the ROI ellipsoid in each dimension.
### Fifth, calculate centromere and background signals
Using these ROIs, average grayscale pixel values are calculated from the data image.
### Sixth, save to file
These data are saved to a .csv file named with the date, time to the minute, and condition. The headers store:
```
Condition: Provided when running the program. Useful when compiling  multiple quantifications.
Cell: The cell number.
Centromere Mean: The average grayscale value in the centromere region.
Centromere Pixels: The number of pixels in the centromere region.
Background Mean: The average grayscale value in the background region.
Background Pixels: The number of pixels in the background region.
Centromere-Background: The difference of Centromere Mean and Background Mean
Normalized: Centromere-Background divided by Background.
```
In most cases, the most reliable data should be in the Centromere-Background column.

ImageJ ROI files are also created to reference back to the original images.

## Tunable Parameters

    centromere_radius: int          # nm; half width at half maximum of centromere signal. Default 150.
    xscale: int                     # nm/pixel; Default 147 for a 63x objective on Leica 2.
    yscale: int                     # nm/pixel; Default 147 for a 63x objective on Leica 2.
    zscale: int                     # nm/slice; Default 500.
    min_centromere_dist: int        # voxels; Chebyshev distance separating centromere center points. Default 2.
    max_num_maxima: int             # -; Maximum number of centromere centers to identify, starting from most intense. Default 38.
    maxima_std_threshold: int       # -; Number of standard deviations above mask image mean intensity to accept centromere point. Default 0.
                                    # NB: including this value enables a significant speed-up, so it is not optional.
    roi_xradius: int or float       # pixels; Radius of centromere ellipsoid in x. Default 4.
    roi_yradius: int or float       # pixels; Radius of centromere ellipsoid in y. Default 4.
    roi_zradius: int or float       # slices; Radius of centromere ellipsoid in z. Default 3.
    self.dist_from_edge_xy: int     # pixels; minimum distance from edge of image for an ROI. Default 20.
    self.dist_from_edge_z: int      # slices; minimum distance from edge of image for an ROI. Default 2.
    bg_xradius: int or float        # pixels; Radius of background ellipsoid in x. Default roi_xradius + 1.
    bg_yradius: int or float        # pixels; Radius of background ellipsoid in y. Default roi_yradius + 1.
    bg_zradius: int or float        # slices; Radius of background ellipsoid in z. Default roi_zradius + 1.
    bg_std_threshold: None or int   # -; Number of standard deviations above quantified background to include data. Ignored if None. Default None.

## Works Cited
Vermolen BJ, Garini Y, Young IT, Dirks RW, Raz V. 2008. "Segmentation and analysis of the three-dimensional 
redistribution of nuclear components in human mesenchymal stem cells." Cytometry 73A: 816-824.
https://doi.org/10.1002/cyto.a.20612.