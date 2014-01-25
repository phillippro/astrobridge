astrobridge
===========

a software to calibrate astro-photos from ccd and cmos cameras

This version can only calibrate JPG images. 

aim: To introduce how to use python to do astronomical image calibration.

function: To calibrate jpg images and generate calibrated jpg and fits images

Required software: alipy (and its required package), numpy, pyraf, scipy, pyfits, Python Imaging Library (PIL) 

1. Make directory

   The main folder: put all stuff, e.g. calibration
   
   The data folder: put photos, e.g. calibration/data
   
   The calibrated folder: put calibrated photos, e.g. calibration/calibrated_image

2. Take photos

   Because our code is aiming at calibrating jpg pictures, please take JPG
astro photos. (Note: Our software can only recognize .JPG rather than .jpg)

   Please put your bias, darks, flats and lights into data/bias, data/darks, data/flats,
data/objects, respectively.


3. Run code

   The bash code is used to calibrate all images.
   You can use it like this: 
   
   source astrobridge_jpg.sh e1 e2
   
   or
   /bin/bash astrobridge_jpg.sh e1 e2
   
   or
   ./astrobridge_jpg.sh e1 e2
   
   , where e1, e2 are exposure time of darks and flats separately. 

   Note: The code astrobridge_jpg.sh and calibration_jpg.sh should be put into calibration/ 


4. stacked image
   
   The stacked image is created using:

   python calibrate_jpg.py stack_image calibrated_image "reference_image_without_extension" "name_of_outfolder"
      
   where the calibrated_image is the folder where the calibrated fits locate. 

   The stacked rgb image is created using: 
   
   python calibrate_jpg.py stack_image create_rgb_image mod
      
   where mod should be "intersect" or "Union". Intersect means the stacked rgb is the intersection between r, g, b fits while the union mode means the sstacked rgb is the union of these files. 


4. Notes: 

   a. Our software will skip the calibratoin of those photos who have already been calibrated (ie the calibrated file exist!)

   b. Our software can be easily updated to calibrate images in .NEF, other raw formats and fits directly.
   
   c. The stacked rgb image may show the displacement between different color frame which may caused by the low astronometric accuracy arising from low quality images (e.g. defocus or sources with star traces). 
   

