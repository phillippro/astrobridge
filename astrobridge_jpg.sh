#!/bin/bash
###############calibrate JPG files
########calibrate_jpg.py is a revised JPG version of Lia's calibration.py
########generate master bias
python calibrate_jpg.py cal_bias data/bias bias_master
########generate master dark
###unscaled dark
#python calibrate_jpg.py cal_dark data/dark dark_master
###scaled dark
python calibrate_jpg.py cal_scale_dark data/dark bias_master scale_dark_master $1 $2
##here e1, e2 are exposure time of darks and flats (lights) separately
##note: In our case, the flat is with an exposure time of 0.01s and the dark is with an exposure time of 5s
########generate master flat fits file ( because JPG cannot store the original values)
python calibrate_jpg.py cal_flat data/flat scale_dark_master flat_master
########calibrate images
DIRECTORY = "calibrated_image"
if [ ! -d $DIRECTORY ]; then
    mkdir calibrated_image
    echo "No folder named calibrated_image!"
fi
for i in $( ls data/object/*.JPG); do
    python calibrate_jpg.py cal_image $i scale_dark_master flat_master "${i%.*}"
done
mv data/object/*.fits data/object/*-calibrate.JPG calibrated_image
#########stack images
###note: here $3 is the reference image, should be input in the command line
###in the test images, I choose ref=DSC0019, i.e. outfolder= stacked_image
if [ $3='' ]; then
    echo 'no reference image!'
else
    python stack_image.py stack_fits calibrated_image $3 "stacked_image"
fi
###combine images to create rgb images, the input parameter is "intersect" or "union"
if [ $4='' ]; then
    echo 'No mode for rgb stack, I will use intersect by default!'
else
    python stack_image.py create_rgb_stacked $4
fi
