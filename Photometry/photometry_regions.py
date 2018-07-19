from astropy.io import fits    #handle fits files
from astropy import wcs    #world coordinate system transformations
from photutils import CircularAperture
from photutils import EllipticalAperture
from photutils import aperture_photometry
import numpy as np
import xlrd
import math

#convert RA/declination from "time"/"arctime" to degrees
def positionStringtoInt(input_string, ra):
    str_split = input_string.split(':') #split the string at the colon
    a, b, c = float(str_split[0]), float(str_split[1]), float(str_split[2])
    int_split = [a, b, c]

    if ra: #perform RA conversion if RA is true
        return(int_split[0] * 15 + int_split[1] * 0.25 + int_split[2] / 240)
    elif int_split[0] > 0: #otherwise, perform declination conversion
        return(int_split[0] + int_split[1] / 60 + int_split[2] / 3600)
    else:
        return(int_split[0] - int_split[1] / 60 - int_split[2] / 3600)

#convert arcsec to number of pixels
def arcsecToPix(arcsec, scale):
    return arcsec / 3600 / scale

#get the solid angle in sr per pixel
def solidAngle(scale):
    return scale * scale * math.pi * math.pi / 180 / 180

#parse the relevant information from a region file line. Return None if the line is not an object.
#Otherwise, return relevant region information in a tuple with a boolean.
#The boolean is True if the aperture is for the galaxy itself, False otherwise
def translate(line, coord_info, scale):
    galaxy_aperture = False
    line_array = line.split('(') #split line at the '(' to determine region shape
    if line_array[0] == 'circle' or line_array[0] == 'ellipse':
        remainder = line_array[1]
        remainder_colorcheck = remainder.split('#') #check if there is a comment on the line
        if len(remainder_colorcheck) == 2:
            if remainder_colorcheck[1] == ' color=red\n': #check if the region is red
                galaxy_aperture = True
        remainder = remainder_colorcheck[0] #now the remaining string looks something like '11:53:59.065,+60:41:01.37,89.7417")'
        remainder_commasplit = remainder.split(',')
        #the first entry of remainder_commasplit is right ascension. Pass directly to the converter.
        ra = positionStringtoInt(remainder_commasplit[0], True)
        #the second entry of remainder_commasplit is declination. Pass directly to the converter.
        dec = positionStringtoInt(remainder_commasplit[1], False)
        #the third entry is the radius for circular apertures, or the semiminor axis for elliptical apertures
        semiminor = float(remainder_commasplit[2].split('"')[0]) #cut the quotation mark away from the number

        #convert ra, dec, and semiminor to pixels
        x, y = coord_info.all_world2pix(ra, dec, 0)
        semiminor_pix = arcsecToPix(semiminor, scale)

        if line_array[0] == 'circle': #if we have a circular object, we're all done!
            return (galaxy_aperture, CircularAperture((x, y), semiminor_pix))
        else: #otherwise, we still need to grab some information
            #the fourth entry is the semimajor axis for elliptical apertures
            semimajor = float(remainder_commasplit[3].split('"')[0]) #cut the quotation mark away from the number
            semimajor_pix = arcsecToPix(semimajor, scale)
            #the fifth entry is the rotation angle for elliptical apertures
            angle = float(remainder_commasplit[4].split(')')[0]) * math.pi / 180 + math.pi / 2 #cut the closing parenthesis away from the number and convert to radians
            return (galaxy_aperture, EllipticalAperture((x, y), semimajor_pix, semiminor_pix, angle))
    else:
        return None

#given the data and a list of background apertures, calculate the mean sky per pixel
def background(value_data, error_data, apertures):
    results = []
    errors = []
    for aperture in apertures: #make a list of flux / pixel and a list of error / pixel
        phot_table = aperture_photometry(value_data, aperture, error_data).as_array()[0] #get photometry results as an array
        results.append(phot_table[3] / aperture.area())
        if len(phot_table) == 4:
            errors.append(0)
        else:
            errors.append(phot_table[4] / aperture.area())
    sum = 0
    for result in results:
        sum += result
    mean = sum / len(results) #calculate mean

    sum = 0
    for error in errors:
        sum += error * error
    mean_error = math.sqrt(sum) / len(errors) #calculate error in mean

    return (mean, mean_error)

#get file paths and load the files
print('\nEnter the file path of the FITS image.')
values_path = str(input('Image File Path: '))

print('\n')

while True:
    try:
        fits_values = fits.open(values_path)
        fits_values.info()
        break
    except ValueError:
        values_path = str(input('Invalid input. Try again: '))
    except FileNotFoundError:
        values_path = str(input('File not found. Try again: '))
    except OSError:
        values_path = str(input('Invalid file. Try again: '))

#print('\n')
#
#print('Enter the level of the FITS file that you would like to analyze. (0, 1, 2,...)')
#while True:
#    try:
#        level = int(input('Level: '))
#        #load image headers
#        values_hdr = fits_values[level].header
#        break
#    except ValueError:
#        print('Invalid input. Try again.')
#    except IndexError:
#        print('This level was not found. Try again.')

level = 1 #for out sample, this is always true
values_hdr = fits_values[level].header

#print('\n')
#
#print('Enter the name of the pixel scale variable in the FITS file. It is usually \'CDELT1\'.')
#while True:
#    try:
#        scale_name = str(input('Variable name: '))
#        pix_scale = abs(values_hdr[scale_name])
#        break
#    except KeyError:
#        print('Variable not found. Try again.')

scale_name = 'CDELT1' # for our sample, this is always true
pix_scale = abs(values_hdr[scale_name])

print('\n')

#parse WCS info from headers
values_wcs = wcs.WCS(values_hdr)

print('Enter the file path of the FITS image error. Enter \'s\' to skip using an error file.')
errors_path = str(input('Error File Path: '))
while True:
    try:
        if errors_path == 's':
            fits_errors = None
            break
        else:
            fits_errors = fits.open(errors_path)
            break
    except ValueError:
        errors_path = str(input('Invalid input. Try again: '))
    except FileNotFoundError:
        errors_path = str(input('File not found. Try again: '))
    except OSError:
        errors_path = str(input('Invalid file. Try again: '))

print('\n')

#retrieve data from the fits files
value_data = fits_values[level].data

if errors_path == 's':
    error_data = None
else:
    error_data = fits_errors[level].data

print('Enter the file path of the regions file. Background apertures should be green, and the object aperture should be red.')
reg_path = str(input('Regions File Path: '))
while True:
    try:
        regions = open(reg_path)
        break
    except ValueError:
        reg_path = str(input('Invalid input. Try again: '))
    except FileNotFoundError:
        reg_path = str(input('File not found. Try again: '))
    except OSError:
        reg_path = str(input('Invalid file. Try again: '))

print('\n')

print('Is this SPIRE data?') #adjust units of SPIRE data. Assumes the fits file has units of MJy / sr for SPIRE and Jy / pix for PACS
spire = ' '
while (spire != 'y' and spire != 'n'):
    spire = str(input('y/n: '))
multiplier = 1
if spire == 'y':
    multiplier = solidAngle(pix_scale) * 1000000 #get the solid angle per pixel, multiply by a million to get Jy

print('\n')

lines = regions.readlines()
background_apertures = []
galaxy_aperture = None
for line in lines: #make the apertures from the regions file and sort into a list of background apertures and a galaxy aperture
    translation = translate(line, values_wcs, pix_scale)
    if translation != None:
        if translation[0]:
            galaxy_aperture = translation[1]
        else:
            background_apertures.append(translation[1])

galaxy_photometry = aperture_photometry(value_data, galaxy_aperture, error_data).as_array()[0] #get photometry results of the galaxy as an array
if len(galaxy_photometry) == 4:
    galaxy_flux = (galaxy_photometry[3], 0) #if no error was calculated, return 0 for error
else:
    galaxy_flux = (galaxy_photometry[3], galaxy_photometry[4])

if len(background_apertures) == 0:
    sky = (0, 0)
else:
    sky = background(value_data, error_data, background_apertures)

final_flux = (galaxy_flux[0] - sky[0] * galaxy_aperture.area()) * multiplier #scale background and subtract
error_sky = (math.sqrt(galaxy_flux[1] * galaxy_flux[1] + sky[1] * galaxy_aperture.area() * sky[1] * galaxy_aperture.area())) * multiplier #calculate error from sky
if spire == 'y': #calculate calibration error
    error_cal = final_flux * 0.07
else:
    error_cal = final_flux * 0.05
final_error = math.sqrt(error_cal * error_cal + error_sky * error_sky)

print('Flux: ' + str(final_flux) + ' Jy\nSky Error: ' + str(error_sky) + ' Jy\nTotal Error: ' + str(final_error) + ' Jy')
