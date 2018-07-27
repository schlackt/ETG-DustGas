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
        #the third entry is the radius for circular apertures, or the one of the axes for elliptical apertures
        radius1 = float(remainder_commasplit[2].split('"')[0]) #cut the quotation mark away from the number

        #convert ra, dec, and radius to pixels
        x, y = coord_info.all_world2pix(ra, dec, 0)
        radius1_pix = arcsecToPix(radius1, scale)

        if line_array[0] == 'circle': #if we have a circular object, we're all done!
            return (galaxy_aperture, CircularAperture((x, y), radius1_pix))
        else: #otherwise, we still need to grab some information
            #the fourth entry is the other axis for elliptical apertures
            radius2 = float(remainder_commasplit[3].split('"')[0]) #cut the quotation mark away from the number
            radius2_pix = arcsecToPix(radius2, scale)
            #the fifth entry is the rotation angle for elliptical apertures
            angle = float(remainder_commasplit[4].split(')')[0]) * math.pi / 180 #cut the closing parenthesis away from the number and convert to radians
            #determine which axis is which
            #the angle reported by ds9 is the angle of radius2 east of north. it needs to be converted to angle of th semimajor axis above the positive x axis
            if radius1 < radius2:
                semiminor = radius1_pix
                semimajor = radius2_pix
                angle = angle + math.pi / 2
            else:
                semiminor = radius2_pix
                semimajor = radius1_pix
            return (galaxy_aperture, EllipticalAperture((x, y), semimajor, semiminor, angle))
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

print('\n')

print('Enter the level of the FITS file that you would like to analyze. (0, 1, 2,...)')
while True:
    try:
        level = int(input('Level: '))
        #load image headers
        values_hdr = fits_values[level].header
        break
    except ValueError:
        print('Invalid input. Try again.')
    except IndexError:
        print('This level was not found. Try again.')

print('\n')

pix_scale = abs(values_hdr['CDELT1']) # degrees / pixel
multiplier = solidAngle(pix_scale) #get the solid angle in steradians per pixel
bmin = abs(values_hdr['BMIN']) # beam minor axis in degrees
bmaj = abs(values_hdr['BMAJ']) # beam major axis in degrees
barea = bmin * bmaj * math.pi * math.pi / 180 / 180 # beam area in sr
bpix = multiplier / barea # beam per pixel

#parse WCS info from headers

# Create a new WCS object.  The number of axes must be set
# from the start
w = wcs.WCS(naxis=2)

# Set up an "Airy's zenithal" projection
# Vector properties may be set with Python lists, or Numpy arrays
w.wcs.crpix = [values_hdr['CRPIX1'], values_hdr['CRPIX2']]
w.wcs.cdelt = np.array([values_hdr['CDELT1'], values_hdr['CDELT2']])
w.wcs.crval = [values_hdr['CRVAL1'], values_hdr['CRVAL2']]
w.wcs.ctype = ["RA---NCP", "DEC--NCP"]

#retrieve data from the fits files
value_data = fits_values[level].data

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

print('Enter the distance to the galaxy in Mpc.')
while True:
    try:
        distance = float(input('Distance (Mpc): '))
        break
    except ValueError:
        distance = float(input('Invalid input. Try again: '))

print('\n')

print('Enter the error in distance to the galaxy in Mpc.')
while True:
    try:
        derror = float(input('Distance Error (Mpc): '))
        break
    except ValueError:
        derror = float(input('Invalid input. Try again: '))

print('\n')

lines = regions.readlines()
background_apertures = []
galaxy_apertures = []
for line in lines: #make the apertures from the regions file and sort into a list of background apertures and a galaxy aperture
    translation = translate(line, w, pix_scale)
    if translation != None:
        if translation[0]:
            galaxy_apertures.append(translation[1])
        else:
            background_apertures.append(translation[1])

#sort the list of galaxy apertures from highest area to lowest area
aperture_areas = []
for aperture in galaxy_apertures:
    aperture_areas.append((aperture.area(), galaxy_apertures.index(aperture)))
aperture_areas_sorted = sorted(aperture_areas, key=lambda aperture: aperture[0]) #sort by area
galaxy_apertures_sorted = []
for aperture in aperture_areas_sorted:
    galaxy_apertures_sorted.insert(0, galaxy_apertures[aperture[1]])

error_data = None

#get background
if len(background_apertures) == 0:
    sky = (0, 0)
else:
    sky = background(value_data, error_data, background_apertures)

#make list of photometry results for each aperture in the form (flux, sky error, total error)
results = []
for galaxy_aperture in galaxy_apertures_sorted:
    galaxy_photometry = aperture_photometry(value_data, galaxy_aperture, error_data).as_array()[0] #get photometry results of the galaxy as an array
    if len(galaxy_photometry) == 4:
        galaxy_flux = (galaxy_photometry[3], 0) #if no error was calculated, return 0 for error
    else:
        galaxy_flux = (galaxy_photometry[3], galaxy_photometry[4])

    final_flux = (galaxy_flux[0] - sky[0] * galaxy_aperture.area()) * bpix #scale background and subtract
    mass = 236000 * distance * distance * final_flux
    error = math.sqrt((2 * 236000 * distance * final_flux * derror)**2 + (mass * 0.15)**2)
    results.append((mass, error))

#now do region subtraction and print results
i = 0
while i < len(results):
    if i == 0:
        print('Galaxy:\n\tHI Mass: ' + str(results[i][0]) + '\n\tError: ' + str(results[i][1]))
    else:
        region_mass = results[i - 1][0] - results[i][0]
        region_error = math.sqrt(results[i - 1][1] * results[i - 1][1] + results[i][1] * results[i][1])
        print('Next Region:\n\tHI Mass: ' + str(region_mass) + '\n\tError: ' + str(region_error))
    if i == len(results) - 1:
        print('Center:\n\tHI Mass: ' + str(results[i][0]) + '\n\tError: ' + str(results[i][1]))
    i += 1
