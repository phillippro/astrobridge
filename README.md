astrobridge
===========

a software to calibrate astro-photos from ccd and cmos cameras

This is version 1.0 which only calibrate JPG images. 

aim: To introduce how to use python to do astronomical image calibration.

function: To calibrate jpg images and generate calibrated jpg and fits images

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


4. Notes: 

   a. Our software will skip the calibratoin of those photos who have already been calibrated (ie the calibrated file exist!)

   b. Our software can be easily updated to calibrate images in .NEF, other raw formats and fits directly.

