from astropy import units as u
from astropy import constants
from dust_emissivity import blackbody, fit_sed
import numpy as np
import matplotlib.pyplot as plt
import csv

_m_p = constants.m_p.cgs.value
_c = 300000000

print('\nEnter the file path of the results spreadsheet.')
results_path = str(input('Spreadsheet File Path: '))

print('\n')

print('Generate plots?')
plots = ' '
while (plots != 'y' and plots != 'n'):
    plots = str(input('y/n: '))

while True:
    try:
        with open(results_path, newline = '') as results:
            reader = csv.reader(results)
            next(reader) # skip header
            data = [row for row in reader]
        break
    except ValueError:
        results_path = str(input('Invalid input. Try again: '))
    except FileNotFoundError:
        results_path = str(input('File not found. Try again: '))
    except OSError:
        results_path = str(input('Invalid file. Try again: '))

print('\n')

writetofile = [] # array that contains the lines of the output csv
for row in data:
    i = 3
    detections = 0
    wavelengths = []
    fluxes = []
    fluxes_err = []
    outputrow = [row[0], row[1], row[2]] # array that contains the elements of a line of the output csv
    while i < 24:
        if row[i] == '': # skip a wavelength with no data
            i += 3
        elif float(row[i]) <= 3 * float(row[i + 1]): # skip a nondetection
            i += 3
        else:
            detections += 1

            if i == 3: # append the wavelength of the detection
                wavelengths.append(24)
            elif i == 6:
                wavelengths.append(70)
            elif i == 9:
                wavelengths.append(100)
            elif i == 12:
                wavelengths.append(160)
            elif i == 15:
                wavelengths.append(250)
            elif i == 18:
                wavelengths.append(350)
            elif i == 21:
                wavelengths.append(500)

            fluxes.append(float(row[i]))
            fluxes_err.append(float(row[i + 2]))

            i += 3

    if detections >= 3:
        maxfreq = _c / (wavelengths[0] * 10**-6)
        minfreq = _c / (wavelengths[len(wavelengths) - 1] * 10**-6)
        wavelengths = np.array(wavelengths) * u.um
        frequencies = wavelengths.to(u.Hz, u.spectral())
        fluxes = np.array(fluxes) * u.Jy
        fluxes_err = np.array(fluxes_err) * u.Jy
        tguess, bguess, nguess, sguess = float(row[24]) * u.K, float(row[25]), float(row[26]) * u.cm**-2, float(row[27])

        t = np.arange(minfreq, maxfreq, minfreq * 0.01)

        try:
            pars, errs, chi = fit_sed.fit_modified_bb(frequencies, fluxes, fluxes_err, (tguess, bguess, nguess, sguess), detections, return_error = True)

            if plots == 'y':
                plt.title(row[0] + ': ' + row[2])
                plt.xlabel('Frequency (Hz)')
                plt.ylabel('Flux (Jy)')
                #plt.plot(frequencies / u.Hz, (fluxes / u.Jy), 'bo', t, (blackbody.modified_blackbody(t * u.Hz, pars[0] * u.K, pars[1], pars[2] * u.cm**-2, pars[3]) / u.erg/u.s/u.cm**2/u.Hz/u.sr) * 10**23, 'r-')
                plt.errorbar(frequencies / u.Hz, (fluxes / u.Jy), yerr = fluxes_err / u.Jy, fmt='ko')
                plt.plot(t, (blackbody.modified_blackbody(t * u.Hz, pars[0] * u.K, pars[1], pars[2] * u.cm**-2, pars[3]) / u.erg/u.s/u.cm**2/u.Hz/u.sr) * 10**23, 'r-')
                plt.show()

                outputrow.append(pars[0])
                outputrow.append(errs[0])
                outputrow.append(pars[1])
                outputrow.append(errs[1])
                outputrow.append(pars[2])
                outputrow.append(errs[2])
                outputrow.append(pars[3])
                outputrow.append(errs[3])
                outputrow.append(chi)
                outputrow.append(detections)
        except ValueError:
            print('Could not generate fit for ' + row[0] + ': ' + row[2] + '. Try changing the guesses.')

    writetofile.append(outputrow)

with open('fitting_results.csv', 'w', newline = '') as fitting_results:
    fitwriter = csv.writer(fitting_results)
    fitwriter.writerow(['Galaxy', 'Aperture', 'ApName', 'Temp', 'Temp Error', 'Beta', 'Beta Error', 'N', 'N Error', 'Scale', 'Scale Error', 'Chi Square','Detection'])
    for row in writetofile:
        fitwriter.writerow(row)


#wavelengths = np.array([24, 100, 160, 250, 350, 500]) * u.um
#frequencies = wavelengths.to(u.Hz, u.spectral())
#flux = np.array([0.18, 0.76, 1.04, 1.74, 0.81, 0.32]) * u.Jy
#flux_error =  np.array([0.03, 0.09, 0.07, 0.12, 0.06, 0.03]) * u.Jy

#tguess, bguess, nguess = 20.*u.K,2.,1e22*u.cm**-2
#bbunit = u.Jy
#pars = fit_sed.fit_modified_bb(frequencies, flux, flux_error, guesses=(tguess, bguess, nguess), return_error = True)
#print(pars)
