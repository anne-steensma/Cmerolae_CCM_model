
## -- load in data set and relevant packages -- ##
install.packages("drc")
install.packages("readxl")
library(dplyr)
library(drc)
library(readxl)

data = read_excel("~/Desktop/co2response_1percent.xlsx")
#put the temperature (in Celsius) that your data was collected at here: v
tem_c = 40

## -- determine offset needed for each rep --
#an "offset" is applied to the data before fitting, 
#because the fitting function does not work if there are negative values in the data.
#later the offset is removed to get the "real" Km

x1 = as.vector(data$rep1co2)
x1 = sort(x1, decreasing = FALSE)
if (x1[1] <= 0) {
  x_offset_1 = abs(x1[1])
} else {
  x_offset_1 = -abs(x1[1]) + -abs(x1[1])*1e-18 #apparently it doesn't work if the point is exactly at (0,0), hence the extremely small addition
}
y1 = as.vector(data$rep1flux)
y1 = sort(y1, decreasing = FALSE)
if (y1[1] <= 0) {
  y_offset_1 = abs(y1[1])
} else {
  y_offset_1 = -abs(y1[1]) + -abs(y1[1])*1e-18 
}

x2 = as.vector(data$rep2co2)
x2 = sort(x2, decreasing = FALSE)
if (x2[1] <= 0) {
  x_offset_2 = abs(x2[1])
} else {
  x_offset_2 = -abs(x2[1]) + -abs(x2[1])*1e-18 #apparently it doesn't work if the point is exactly at (0,0), hence the extremely small addition
}
y2 = as.vector(data$rep2flux)
y2 = sort(y2, decreasing = FALSE)
if (y2[1] <= 0) {
  y_offset_2 = abs(y2[1])
} else {
  y_offset_2 = -abs(y2[1]) + -abs(y2[1])*1e-18 
}

x3 = as.vector(data$rep3co2)
x3 = sort(x3, decreasing = FALSE)
if (x3[1] <= 0) {
  x_offset_3 = abs(x3[1])
} else {
  x_offset_3 = -abs(x3[1]) + -abs(x3[1])*1e-18 #apparently it doesn't work if the point is exactly at (0,0), hence the extremely small addition
}
y3 = as.vector(data$rep3flux)
y3 = sort(y3, decreasing = FALSE)
if (y3[1] <= 0) {
  y_offset_3 = abs(y3[1])
} else {
  y_offset_3 = -abs(y3[1]) + -abs(y3[1])*1e-18 
}



## -- apply offsets to the data --

data_offset = mutate(data, rep1co2 = rep1co2 + x_offset_1)
data_offset = mutate(data_offset, rep2co2 = rep2co2 + x_offset_2)
data_offset = mutate(data_offset, rep3co2 = rep3co2 + x_offset_3)

data_offset = mutate(data_offset, rep1flux = rep1flux + y_offset_1)
data_offset = mutate(data_offset, rep2flux = rep2flux + y_offset_2)
data_offset = mutate(data_offset, rep3flux = rep3flux + y_offset_3)

# -- perform Michaelis-Menten fitting on offset data from each rep --

S1 = as.vector(data_offset$rep1co2)
v1 = as.vector(data_offset$rep1flux)
mrep1 = drm(v1 ~ S1, data = data_offset, fct = MM.2())

S2 = as.vector(data_offset$rep2co2)
v2 = as.vector(data_offset$rep2flux)
mrep2 = drm(v2 ~ S2, data = data_offset, fct = MM.2())

S3 = as.vector(data_offset$rep3co2)
v3 = as.vector(data_offset$rep3flux)
mrep3 = drm(v3 ~ S3, data = data_offset, fct = MM.2())

# -- store model parameters and summary (including Km in terms of estimated dissolved CO2) --

Km_offset = c(mrep1$coefficients[2], mrep2$coefficients[2], mrep3$coefficients[2])
Vmax_offset = c(mrep1$coefficients[1], mrep2$coefficients[1], mrep3$coefficients[1])

offset_x = c(x_offset_1, x_offset_2, x_offset_3)
offset_y = c(y_offset_1, y_offset_2, y_offset_3)

Km_real = c(mrep1$coefficients[2] - x_offset_1, mrep2$coefficients[2] - x_offset_2, mrep3$coefficients[2] - x_offset_3)
Vmax_real = c(mrep1$coefficients[1] - y_offset_1, mrep2$coefficients[1] - y_offset_2, mrep3$coefficients[1] - y_offset_3)

#(function that calculates liquid-phase Km from gas-phase Km):

dissolved_Km = function(x) {
  tem = tem_c+273.15 #parameter: temperature in K
  tem_c = tem - 273.15 #PARAMETER: temperature in C
  ppm_CO2 = x  #PARAMETER: atmospheric CO2 concentration in ppm
  air_pressure = 101325 #PARAMETER: air pressure in Pascal (hPa * 100)
  pp_CO2 = (ppm_CO2/1000000) * air_pressure * (1/100000) #calculates: partial pressure of CO2 in bar (see pv = nrt, or p = nrt)
  
  henry_constant_std = 0.035 #PARAMETER: henry's law constant in mol/(kg*bar)
  tem_dependence_constant = 2400 #PARAMETER: henry's law temperature dependence constant in K
  henry_constant_new = henry_constant_std * 2.718^((tem_dependence_constant)*((1/tem) - (1/298.15))) #calculates: temperature-adjusted henry's constant in mol/(kg*bar)
  
  CO2_conc_molkg = pp_CO2 * henry_constant_new #calculates: CO2 concentration in mol/kg
  water_density = (((-5e-06)*(tem_c^2)) + (8e-06*tem_c) + 1.0001)*10e02 #calculates: water density in kg/m^3 at temperature of interest (should be ~1000)
  CO2_conc_molar = water_density * CO2_conc_molkg * 0.001 #calculates: molar concentration of total DIC
}

Km_real_aq = c(dissolved_Km(Km_real[1]), dissolved_Km(Km_real[2]), dissolved_Km(Km_real[3]))

# -- create dataframe --
parameters_df = data.frame(Km_offset, Vmax_offset, offset_x, offset_y, Km_real, Vmax_real, Km_real_aq)


