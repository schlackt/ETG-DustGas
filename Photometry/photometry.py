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
    str_split = input_string.split(' ') #split the string at the spaces

    #exception handling
    while True:
        try:
            while len(str_split) != 3:
                input_string = str(input('Incorrect format. Try again: '))
                str_split = input_string.split(' ')
            a, b, c = float(str_split[0]), float(str_split[1]), float(str_split[2])
            break
        except ValueError:
            input_string = str(input('Invalid number. Try again: '))
            str_split = input_string.split(' ')

    int_split = [a, b, c]

    if ra: #perform RA conversion if RA is true
        return(int_split[0] * 15 + int_split[1] * 0.25 + int_split[2] / 240)
    else: #otherwise, perform declination conversion
        return(int_split[0] + int_split[1] /60 + int_split[2] / 3600)

#convert arcsec to number of pixels
def arcsecToPix(arcsec, scale):
    return arcsec / 3600 / scale

#perform photometry on a circular aperture. Returns tuple in the form (x position, y position, aperture sum, error)
def circularPhotometry(value_data, error_data, coord_info, scale, ra, dec, radius):
    # convert to pixels
    x, y = coord_info.all_world2pix(ra, dec, 0)
    radius_pix = arcsecToPix(radius, scale)

    #perform photometry
    aperture = CircularAperture((x, y), radius_pix)
    phot_table = aperture_photometry(value_data, aperture, error_data)

    results = phot_table.as_array()[0] #pull results from table

    if len(results) == 4:
        return (results[1], results[2], results[3], 0) #if no error was calculated, return 0 for error
    else:
        return (results[1], results[2], results[3], results[4])


def ellipticalPhotometry(value_data, error_data, coord_info, scale, ra, dec, semimajor, semiminor, angle):
    # convert to pixels
    x, y = coord_info.all_world2pix(ra, dec, 0)
    semimajor_pix = arcsecToPix(semimajor, scale)
    semiminor_pix = arcsecToPix(semiminor, scale)

    #perform photometry
    aperture = EllipticalAperture((x, y), semimajor_pix, semiminor_pix, angle)
    phot_table = aperture_photometry(value_data, aperture, error_data)

    results = phot_table.as_array()[0] #pull results from table

    if len(results) == 4:
        return (results[1], results[2], results[3], 0) #if no error was calculated, return 0 for error
    else:
        return (results[1], results[2], results[3], results[4])

#get file paths and load the files
print('Enter the file path of the FITS image.')
values_path = str(input('Image File Path: '))
while True:
    try:
        fits_values = fits.open(values_path)
        break
    except ValueError:
        values_path = str(input('Invalid input. Try again: '))
    except FileNotFoundError:
        values_path = str(input('File not found. Try again: '))
    except OSError:
        values_path = str(input('Invalid file. Try again: '))

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

#get right ascension and declination
ra_in = str(input('Enter right ascension of the object in the format \'hour min sec\': '))
ra = positionStringtoInt(ra_in, True)

dec_in = str(input('Enter declination of the object in the format \'deg arcmin arcsec\': '))
dec = positionStringtoInt(dec_in, False)

#get aperture type
print('Would you like to use a circular or elliptical aperture?')
aperture_type = ' '
while (aperture_type != 'circular' and aperture_type != 'elliptical'):
    aperture_type = str(input('Enter \'circular\' or \'elliptical\': '))

#ask about sky subtraction
print('Would you like to perform a background subtraction?')
sky_subtraction = ' '
while (sky_subtraction != 'y' and sky_subtraction != 'n'):
    sky_subtraction = str(input('y/n: '))



#load image headers
values_hdr = fits_values[1].header
pix_scale = abs(values_hdr['CDELT1'])

#parse WCS info from headers
values_wcs = wcs.WCS(values_hdr)

#retrieve data from the fits files
value_data = fits_values[1].data

if errors_path == 's':
    error_data = None
else:
    error_data = fits_errors[1].data

#get the necessary information for the selected aperture
if aperture_type == 'circular':
    while True:
        try:
            radius = abs(float(input('Enter the aperture radius in arcseconds: ')))
            break
        except ValueError:
            print('Invalid input.')

    photometry_results = circularPhotometry(value_data, error_data, values_wcs, pix_scale, ra, dec, radius)
else:
    while True:
        try:
            semimajor = abs(float(input('Enter the semimajor axis in arcseconds: ')))
            semiminor = abs(float(input('Enter the semiminor axis in arcseconds: ')))
            angle = float(input('Enter the rotation angle in degrees: ')) * math.pi / 180
            break
        except ValueError:
            print('Invalid input.')

    photometry_results = ellipticalPhotometry(value_data, error_data, values_wcs, pix_scale, ra, dec, semimajor, semiminor, angle)

print(photometry_results)
fits_values.close()
