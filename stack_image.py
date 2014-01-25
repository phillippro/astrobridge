#! /usr/bin/env python
import pyfits
import numpy as np
import matplotlib.pyplot as plt
from calibrate_jpg import Image
from scipy.misc import imread, imsave
from astropy.io import fits
from os import listdir
import os, sys
#from os import glob
import glob
import Image as Img
import alipy
##note: alipy depend on pyraf, thus if you cannot install pyraf successfully, you cannot use alipy.align.irafalign, you can only use alipy.align.affineremap to stack. The other simple way to install pyraf is to use Ureka
###################################
#####stack functions
###################################
def stackimage( img_list,ref_image, outfile):
    images_to_align = img_list
    identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
# Put visu=True to get visualizations in form of png files (nice but much slower) On multi-extension data, you will want to specify the hdu (see API doc).
    for id in identifications:
        if id.ok == True: # i.e., if it worked
                print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
                # id.trans is a alipy.star.SimpleTransform object. Instead of printing it out as a string,
                # you can directly access its parameters :
                #print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
                #print id.trans.matrixform()
                #print id.trans.inverse() # this returns a new SimpleTransform object
        else:
                print "%20s : no transformation found !" % (id.ukn.name)

    outputshape = alipy.align.shape(ref_image)
    for id in identifications:
        if id.ok == True:
                # Variant 1, using only scipy and the simple affine transorm :
              #  alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, alifilepath=outfile, makepng=True)

                # Variant 2, using geomap/gregister, correcting also for distortions :
                alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape, alifilepath=outfile, makepng=False)
                # id.uknmatchstars and id.refmatchstars are simply lists of corresponding Star objects.
                # By default, the aligned images are written into a directory "alipy_out".

def crop_image(inlist, outlist, mod):
    ##mod=intersect or union
    ##img_list, the r, g, b stacked files
    img0 = Image(inlist[0])
    img1 = Image(inlist[1])
    img2 = Image(inlist[2])
    if mod=="intersect":
        ind_zero = ((img0.data <=1.0) | (img1.data <=1.0) | (img2.data <=1.0))
    elif mod == "union":
        ind_zero = ((img0.data <=1.0) & (img1.data <=1.0) & (img2.data <=1.0))
    else:
        ind_zero = (img0.data <=1.0)
        
    ind_rows = np.all(ind_zero, axis=0)
    ind_cols = np.all(ind_zero, axis=1)
    x_ind_to_rm = np.where(ind_cols==True)[0]
    y_ind_to_rm = np.where(ind_rows==True)[0]
    
    tmp = img0.data[:,[y for y in range(img0.data.shape[1]) if y not in y_ind_to_rm]]
    img0.data = tmp[[x for x in range(img0.data.shape[0]) if x not in x_ind_to_rm],:]
    img0.write_to(outlist[0])
    tmp = img1.data[:,[y for y in range(img1.data.shape[1]) if y not in y_ind_to_rm]]
    img1.data = tmp[[x for x in range(img1.data.shape[0]) if x not in x_ind_to_rm],:]
    img1.write_to(outlist[1])
    tmp = img2.data[:,[y for y in range(img2.data.shape[1]) if y not in y_ind_to_rm]]
    img2.data = tmp[[x for x in range(img2.data.shape[0]) if x not in x_ind_to_rm],:]
    img2.write_to(outlist[2])
    return

###########################################
###executable function without parameters
###########################################
def stack_fits():
    # inputs: calibrated images, reference image
    imagefolder, ref, outfolder = sys.argv[2:5]
    for i in range(0,3):
        if i==0:
            imagefile  = [ FF for FF in glob.glob(imagefolder+'/*-r.fits') ]
            outfile = outfolder+'/stacked-r.fits'
            ref_image = imagefolder+'/'+ref+'-r.fits'
        elif i==1:
            imagefile  = [ FF for FF in glob.glob(imagefolder+'/*-g.fits') ]
            outfile = outfolder+'/stacked-g.fits'
            ref_image = imagefolder+'/'+ref+'-r.fits'
        else:
            imagefile  = [ FF for FF in glob.glob(imagefolder+'/*-b.fits') ]
            outfile = outfolder+'/stacked-b.fits'
            ref_image = imagefolder+'/'+ref+'-r.fits'
            
# Put visu=True to get visualizations in form of png files (nice but much slower) On multi-extension data, you will want to specify the hdu (see API doc).    
        stackimage(  imagefile, ref_image, outfile)

    
def create_rgb_stacked():
#step1: align r, g, b images
    mod = sys.argv[2]
    imgr = 'stacked_image/stacked-r.fits'
    imgg = 'stacked_image/stacked-g.fits'
    imgb = 'stacked_image/stacked-b.fits'
    inlist = [imgr,imgg,imgb]
    outlist = ['stacked_image/stacked0.fits','stacked_image/stacked1.fits','stacked_image/stacked2.fits']#output media frames
    print "mod name:", mod
    crop_image(inlist,outlist,mod)
    img_rgb = outlist
    ref_image = outlist[0]#randomly selected
    outputshape = alipy.align.shape(ref_image)
    print "outputshape:",outputshape
    identifications = alipy.ident.run(ref_image, img_rgb, visu=False)
    out_no = [0,1,2]
    for i in range(0,3):
        id = identifications[i]
        if id.ok == True:
            outfile = "stacked_image/stacked-"+mod+str(out_no[i])+".fits"#output final fits
#            alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, alifilepath=outfile, makepng=True)
            alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape, alifilepath=outfile, makepng=False)    
        
    img0 = Image('stacked_image/stacked-'+mod+'0.fits')
    img1 = Image('stacked_image/stacked-'+mod+'1.fits')
    img2 = Image('stacked_image/stacked-'+mod+'2.fits')
    img_final = Image('')
    img_final.data = np.dstack([img0.data,img1.data,img2.data])
    imsave( 'stacked_image/stacked-'+mod+'-rgb.JPG', img_final.data )
    
# https://gist.github.com/adrn/6690910
if __name__ == "__main__":
    func = getattr(sys.modules[__name__], sys.argv[1])
    func()
