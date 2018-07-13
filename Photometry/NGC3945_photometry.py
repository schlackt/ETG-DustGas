from photutils import CircularAperture
from photutils import aperture_photometry
import numpy as np
import xlrd

position350 = (95, 94)
radius350 = 17

position250 = (155.5, 154)
radius250 = 27.826

position160 = (118.5, 144)
radius160 = 21.205

#create the aperture
aperture350 = CircularAperture(position350, radius350)
aperture250 = CircularAperture(position250, radius250)
aperture160 = CircularAperture(position160, radius160)

#open spreadsheets with image and error data
data_workbook = xlrd.open_workbook('D:\\ETG\\Data\\Photometry\\spreadsheets\\NGC3945.xlsx')
data_worksheet_350 = data_workbook.sheet_by_name('350')
data_worksheet_250 = data_workbook.sheet_by_name('250')
data_worksheet_250_conv = data_workbook.sheet_by_name('250_conv')
data_worksheet_160 = data_workbook.sheet_by_name('160')
data_worksheet_160_conv = data_workbook.sheet_by_name('160_conv')

error_workbook = xlrd.open_workbook('D:\\ETG\\Data\\Photometry\\spreadsheets\\NGC3945_error.xlsx')
error_worksheet_350 = error_workbook.sheet_by_name('350')
error_worksheet_250 = error_workbook.sheet_by_name('250')
error_worksheet_250_conv = error_workbook.sheet_by_name('250_conv')
error_worksheet_160 = error_workbook.sheet_by_name('160')
error_worksheet_160_conv = error_workbook.sheet_by_name('160_conv')

#create empty 2D arrays of the proper size
rows350 = 190
columns350 = 189
data350 = [[0 for x in range(columns350)] for y in range(rows350)]
error350 = [[0 for x in range(columns350)] for y in range(rows350)]

rows250 = 313
columns250 = 308
data250 = [[0 for x in range(columns250)] for y in range(rows250)]
error250 = [[0 for x in range(columns250)] for y in range(rows250)]
data250_conv = [[0 for x in range(columns250)] for y in range(rows250)]
error250_conv = [[0 for x in range(columns250)] for y in range(rows250)]

rows160 = 238
columns160 = 289
data160 = [[0 for x in range(columns160)] for y in range(rows160)]
error160 = [[0 for x in range(columns160)] for y in range(rows160)]
data160_conv = [[0 for x in range(columns160)] for y in range(rows160)]
error160_conv = [[0 for x in range(columns160)] for y in range(rows160)]

#put the data and errors into the 2D arrays for Photometry
for row in range(0, rows350): #loop over the rows
    for column in range(0, columns350): #loop over the columns
        data350[row][column] = data_worksheet_350.cell(row, column).value
        error350[row][column] = error_worksheet_350.cell(row, column).value

for row in range(0, rows250): #loop over the rows
    for column in range(0, columns250): #loop over the columns
        data250[row][column] = data_worksheet_250.cell(row, column).value
        error250[row][column] = error_worksheet_250.cell(row, column).value
        data250_conv[row][column] = data_worksheet_250_conv.cell(row, column).value
        error250_conv[row][column] = error_worksheet_250_conv.cell(row, column).value

for row in range(0, rows160): #loop over the rows
    for column in range(0, columns160): #loop over the columns
        data160[row][column] = data_worksheet_160.cell(row, column).value
        error160[row][column] = error_worksheet_160.cell(row, column).value
        data160_conv[row][column] = data_worksheet_160_conv.cell(row, column).value
        error160_conv[row][column] = error_worksheet_160_conv.cell(row, column).value

#do the photometry
phot_table350 = aperture_photometry(data350, aperture350, error350)
phot_table250 = aperture_photometry(data250, aperture250, error250)
phot_table250_conv = aperture_photometry(data250_conv, aperture250, error250_conv)
phot_table160 = aperture_photometry(data160, aperture160, error160)
phot_table160_conv = aperture_photometry(data160_conv, aperture160, error160_conv)
print("350 um Photometry:")
print(phot_table350)
print(" ")
print("250 um Photometry:")
print(phot_table250_conv)
print(" ")
print("250 um Photometry (unconvolved):")
print(phot_table250)
print(" ")
print("160 um Photometry:")
print(phot_table160_conv)
print(" ")
print("160 um Photometry (unconvolved):")
print(phot_table160)
print(" ")
