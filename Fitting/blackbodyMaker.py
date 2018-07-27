from astropy import units as u
from dust_emissivity import blackbody
import numpy as np
wavelengths = np.array([100,160]) * u.um
frequencies = wavelengths.to(u.Hz, u.spectral())
temperatures = np.array([10,20,30,40]) * u.K
column = 1e22 * u.cm**-2

for temperature in temperatures:
    flux = blackbody.modified_blackbody(frequencies, temperature, beta=2)
    ratio = flux[0] / flux[1]
    print('Temperature: ' + str(temperature))
    print('Ratio: ' + str(ratio))
    print('\n')
