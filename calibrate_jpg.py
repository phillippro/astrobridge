#! /usr/bin/env python
import pyfits
import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import imread, imsave
from astropy.io import fits
from os import listdir
import os, sys
#from os import glob
import glob
import Image as Img
import alipy

#--------------- Supporting objects ---------------#

class Image(object):
    def __init__( self, filename ):
        
        self.filename = filename
        if filename == '': 
            self.data = 0.0
        elif ".JPG" in filename:
            jpg_image_data = imread( filename ).astype(np.float32)
            self.data = jpg_image_data
        elif ".fits" in filename:
            hdu_list = fits.open( filename )
            self.data = hdu_list[0].data
            hdu_list.close()
        else:
             self.data = 0.0

    def write_to( self, outfile ):
        hdu = fits.PrimaryHDU( self.data )
        hdulist = fits.HDUList([hdu])
        hdulist.writeto( outfile,clobber=True )
        return

#--------------- Image tasks ----------------------#

# http://stackoverflow.com/questions/13753251/median-combining-fits-images-in-python
def combine_image( img_list, combine_mode='median' ):
    result_concat = []
    for ffile in img_list:
        print ffile
        try:
            temp = Image(ffile).data
            result_concat.append( temp )
        except:
            pass
    
    ## This works if every header is correct
    #result_concat = np.ndarray( [ Image(file).data for file in img_list ] )
    
    result = Image( '' )
    if combine_mode == 'median':
        result.data = np.median( result_concat, axis=0 )
    elif combine_mode == 'mean':
        result.data = np.mean( result_concat, axis=0 )
    else:
        result.data = np.zeros( shape = (result.xlen, result.ylen) )
    return result

def scaledark( img_list, biasfile, e1, e2, combine_mode='median' ):
    comb_dark = combine_image( img_list, combine_mode)
    bias_red = Image(biasfile+'-r.fits')
    bias_green = Image(biasfile+'-g.fits')
    bias_blue = Image(biasfile+'-b.fits')
    bias = Image('')
    bias.data = np.dstack([bias_red.data,bias_green.data,bias_blue.data])
    clear_dark = comb_dark.data - bias.data
#    print type(e1), type(e2)
    scale_clear_dark = clear_dark*(float(e2)/float(e1))
    scale_dark = scale_clear_dark + bias.data
    result      = Image('')
    result.data = scale_dark
    return result

def combine_flats( img_list, darkfile, **kwargs ):
    flat_list = []
    dark_red = Image(darkfile+'-r.fits')
    dark_green = Image(darkfile+'-g.fits')
    dark_blue = Image(darkfile+'-b.fits')
    dark = Image('')
    dark.data = np.dstack([dark_red.data,dark_green.data,dark_blue.data])
    for i in img_list:
        print i
        flat_single = Image(i)
        flat_clear = flat_single.data - dark.data
        flat_norm = flat_clear/np.mean( flat_clear )
        try:
            flat_list.append( flat_norm )
        except:
            pass

    result = Image( '' )
    result.data = np.median(flat_list,axis=0)
#    print "length of flat list: ", np.len(flat)
#    print result.data
#    print "mean of red flat: ", np.mean(result.data[:,:,0])
#    print "mean of green flat: ", np.mean(result.data[:,:,1])
#    print "mean of blue flat: ", np.mean(result.data[:,:,2])
    return result

def calibrate_image( img_filename, darkfile, flatfile, outfile=None ):
    image  = Image( img_filename )
    dark   = Image( '' )
    flat   = Image( '' )
    dark.data = np.dstack([Image(darkfile+'-r.fits').data,Image(darkfile+'-g.fits').data,Image(darkfile+'-b.fits').data])
    flat.data = np.dstack([Image(flatfile+'-r.fits').data,Image(flatfile+'-g.fits').data,Image(flatfile+'-b.fits').data])
    
    result = Image('')
    result.data = ( image.data - dark.data )
    result.data = result.data / flat.data
#    goodpix = np.where( flat.data >1.e-4)[0]
#    badpix = np.where( flat.data <=1.e-4)[0]
#    result.data[goodpix] = result.data[goodpix] / flat.data[goodpix]
#    result.data[badpix]  = 1.e-4
    return result

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
                alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, alifilepath=outfile, makepng=True)

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

#######return stacked image
#    result = Image('')
#    result.data = ( image.data - dark.data )
#    result.data = result.data / flat.data
    return 

#------------- For calling from bash shell ------------#
def cal_bias():
    # inputs: darkfolder, outfile
    biasfolder, outfile = sys.argv[2:4]
    biaslist  = [ FF for FF in glob.glob(biasfolder+'/*.JPG') ]
    superbias = combine_image( biaslist )
    print "mean red bias: ", np.mean(superbias.data[:,:,0])
    superbias0 = Image('')
    superbias1 = Image('')
    superbias2 = Image('')
    superbias0.data = superbias.data[:,:,0]
    superbias1.data = superbias.data[:,:,1]
    superbias2.data = superbias.data[:,:,2]
    superbias0.write_to( outfile+"-r.fits" )
    superbias1.write_to( outfile+"-g.fits" )
    superbias2.write_to( outfile+"-b.fits" )
    print "test:"
    print "new mean red bias:", np.mean(Image( outfile+"-r.fits" ).data)
    return

def cal_dark():
    # inputs: darkfolder, outfile
    darkfolder, outfile = sys.argv[2:4]
    darklist  = [ FF for FF in glob.glob(darkfolder+'/*.JPG') ]
    superdark = combine_image( darklist )
    superdark0 = Image('')
    superdark1 = Image('')
    superdark2 = Image('')
    superdark0.data = superdark.data[:,:,0]
    superdark1.data = superdark.data[:,:,1]
    superdark2.data = superdark.data[:,:,2]
    superdark0.write_to( outfile+"-r.fits" )
    superdark1.write_to( outfile+"-g.fits" )
    superdark2.write_to( outfile+"-b.fits" )
    return

def cal_scale_dark():
    # inputs: darkfolder, outfile
    darkfolder, biasfile, outfile, e1, e2 = sys.argv[2:7]
    darklist  = [ FF for FF in glob.glob(darkfolder+'/*.JPG') ]
    scale_dark = scaledark( darklist, biasfile, e1, e2)
    scale_dark0 = Image('')
    scale_dark1 = Image('')
    scale_dark2 = Image('')
    scale_dark0.data = scale_dark.data[:,:,0]
    scale_dark1.data = scale_dark.data[:,:,1]
    scale_dark2.data = scale_dark.data[:,:,2]
    scale_dark0.write_to( outfile+"-r.fits" )
    scale_dark1.write_to( outfile+"-g.fits" )
    scale_dark2.write_to( outfile+"-b.fits" )
    return

def cal_flat():
    # inputs: flatfolder, darkfile, outfile
    flatfolder, darkfile, outfile = sys.argv[2:5]
    flatlist  = [ FF for FF in glob.glob(flatfolder+'/*.JPG') ]
    superflat = combine_flats( flatlist, darkfile )
    superflat0 = Image('')
    superflat1 = Image('')
    superflat2 = Image('')
    superflat0.data = superflat.data[:,:,0]
    superflat1.data = superflat.data[:,:,1]
    superflat2.data = superflat.data[:,:,2]
    superflat0.write_to( outfile+"-r.fits" )
    superflat1.write_to( outfile+"-g.fits" )
    superflat2.write_to( outfile+"-b.fits" )
    return

def cal_image():
    # inputs: imagefile, darkfile, flatfile, outfile
    imagefile, darkfile, flatfile, outfile = sys.argv[2:6]
    final_image = calibrate_image( imagefile, darkfile, flatfile )
    final_image0 = Image('')
    final_image1 = Image('')
    final_image2 = Image('')
    final_image0.data = final_image.data[:,:,0]
    final_image1.data = final_image.data[:,:,1]
    final_image2.data = final_image.data[:,:,2]
    final_image0.write_to( outfile+"-r.fits" )
    final_image1.write_to( outfile+"-g.fits" )
    final_image2.write_to( outfile+"-b.fits" )
    imsave( outfile + '-calibrate.JPG', final_image.data )

def stack_image():
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

    
def create_rgb_image():
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

