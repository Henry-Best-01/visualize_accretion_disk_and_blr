import numpy as np
from amoeba.Classes.accretion_disk import AccretionDisk
from amoeba.Classes.blr_streamline import Streamline
from amoeba.Classes.blr import BroadLineRegion
from amoeba.Util.util import (
    create_maps,
    calculate_gravitational_radius,
    accretion_disk_temperature
)
import matplotlib.pyplot as plt
import streamlit as st
import glob
from astropy.io import fits


path_to_raytraces = "data/"

st.title("The Accretion Disk")
st.write("This is a toy GUI designed to explore how the flux distribution of the accretion disk is related to various parameters")




left_col, right_col = st.columns(2)

mexp = left_col.slider("mass exponent", min_value=6.0, max_value=10.0, step=0.1, value=8.0)
redshift = right_col.slider("redshift", min_value=0.0, max_value=9.0, step=0.1, value=1.0)
inclination = left_col.slider("inclination angle [deg.]", min_value=1, max_value=89, step=1, value=20)
edd_ratio = left_col.slider("Eddington ratio", min_value=0.01, max_value=0.3, step=0.01, value=0.10)
wind_beta = right_col.slider("wind strength", min_value=0.0, max_value=0.8, step=0.01, value=0.0)

axis_range = right_col.slider(r"axis range [$r_{\rm{g}}$]", min_value=10, max_value=1000, step=10, value=100)
apply_gr = st.toggle("apply GR")
wavelength = st.slider("observer frame wavelength range [nm]", min_value=150, max_value=2000, step=5, value=(400, 600))

# grab the GR-ray trace
fname = path_to_raytraces+"RayTrace"+str(int(inclination)).zfill(2)+".fits"
with fits.open(fname) as f:
    r_map = f[0].data
    phi_map = f[1].data
    g_map = f[2].data
    header = f[0].header
    
# work on conversion to mags
min_wavelength = np.min(wavelength) * 10**-9
max_wavelength = np.max(wavelength) * 10**-9
min_frequency = 3e8 / max_wavelength
max_frequency = 3e8 / min_wavelength
delta_freq = abs(max_frequency - min_frequency)
delta_lam = abs(max_wavelength - min_wavelength)


# do some amoeba construction
acc_disk_dict = create_maps(
    mexp,
    redshift_source=redshift+0.0001,
    number_grav_radii=header['numgrs'],
    inclination_angle=inclination,
    resolution=np.size(r_map, 0),
    eddington_ratio=edd_ratio,
    visc_temp_prof="SS",
    temp_beta=wind_beta
)

# adjust maps to include GR
if apply_gr:
    acc_disk_dict['radii_array'] = r_map
    acc_disk_dict['phi_array'] = phi_map
    acc_disk_dict['g_array'] = g_map
    grav_rad = calculate_gravitational_radius(10**mexp)
    t_map = accretion_disk_temperature(r_map * grav_rad, 6.0 * grav_rad, 10**mexp, edd_ratio, beta=wind_beta)
    acc_disk_dict['temp_array'] = t_map


disk = AccretionDisk(**acc_disk_dict)

emission_array = disk.calculate_surface_intensity_map(np.mean(wavelength))
flux_array = emission_array.flux_array
total_flux = emission_array.total_flux
flux_exp = round(np.log10(total_flux), 0)
X, Y = emission_array.get_plotting_axes()


lum_dist = disk.lum_dist
approx_ab_mag = round(-2.5 * np.log10(total_flux / (4 * np.pi * lum_dist**2) * delta_lam / delta_freq * 1000) - 48.6, 1)

title_string = "total emitted flux: "+str(total_flux)[:4]+"e"+str(flux_exp)[:-2]+" W/m"+r", AB mag $\approx$"+str(approx_ab_mag)

fig, ax = plt.subplots(figsize=(8, 6))
contours = ax.contourf(X, Y, (flux_array), 50)
cbar = plt.colorbar(contours, ax=ax, label=r'flux density [W/m$^{2}$/m]')
ax.set_xlabel(r"x [$r_{\rm{g}}$]")
ax.set_ylabel(r"y [$r_{\rm{g}}$]")
ax.set_title(title_string)
ax.set_xlim(-axis_range, axis_range)
ax.set_ylim(-axis_range, axis_range)
ax.set_aspect(1)


st.write(fig)

fig2, ax2 = plt.subplots(figsize=(8, 3))
ax2.plot(X[0, :], flux_array[500, :])
ax2.set_xlim(-axis_range, axis_range)
ax2.set_xlabel(r"x [$r_{\rm{g}}$]")
ax2.set_ylabel(r"flux density across center [W/m$^{2}$/m]")

st.write(fig2)





st.title("The BLR")

left_col, right_col = st.columns(2)

blr_radii = left_col.slider(r"wind radius [$r_{\rm{g}}$]", min_value=10, max_value=2000, step=10, value=(100, 500))
blr_angles = right_col.slider("launch angles [deg.]", min_value=1, max_value=80, step=1, value=(20, 45))
blr_characteristic_distance_1 = left_col.slider(r"inner characteristic distance [$r_{\rm{g}}$]", min_value=10, max_value=2000, step=10, value=200)
blr_characteristic_distance_2 = right_col.slider(r"outer characteristic distance [$r_{\rm{g}}$]", min_value=10, max_value=2000, step=10, value=200)
asympt_vel_1 = left_col.slider("inner asymptotic velocity [c]", min_value=0.001, max_value=0.4, step=0.001, value=0.05)
asympt_vel_2 = right_col.slider("outer asymptotic velocity [c]", min_value=0.001, max_value=0.4, step=0.001, value=0.05)

blr_opt_radius = left_col.slider(r"optimal blr emission radius [$r_{\rm{g}}$]", min_value=10, max_value=4000, step=10, value=300)
blr_opt_width = right_col.slider(r"optimal blr emission width [$r_{\rm{g}}$]", min_value=10, max_value=1000, step=10, value=50)
emitted_wavelength = st.slider("rest frame emission line [nm]", min_value=50, max_value=2000, step=1, value=210)


# Do the Amoeba production of a BLR object
stream_1 = Streamline(blr_radii[0], blr_angles[0], 1000, blr_characteristic_distance_1+0.001, asympt_vel_1, height_step=20)
stream_2 = Streamline(blr_radii[1], blr_angles[1], 1000, blr_characteristic_distance_2+0.001, asympt_vel_2, height_step=20)

my_blr = BroadLineRegion(mexp, 1000, emitted_wavelength, redshift, radial_step=20, height_step=20)
my_blr.add_streamline_bounded_region(stream_1, stream_2)

# define some emission efficiency array based on the optimal emission radius
R, Z = my_blr.get_density_axis()
r_spherical = np.sqrt(R**2 + Z**2)
efficiency = np.exp(-(r_spherical-blr_opt_radius)**2/(2 * blr_opt_width**2))

blr_projection = my_blr.project_blr_intensity_over_velocity_range(inclination, observed_wavelength_range_in_nm=wavelength, emission_efficiency_array=efficiency)

X, Y = blr_projection.get_plotting_axes()
blr_flux = blr_projection.flux_array

fig3, ax3 = plt.subplots(figsize=(8, 6))
contours3 = ax3.contourf(X, Y, blr_flux, 50)
cbar3 = plt.colorbar(contours3, ax=ax3, label="relative blr emission")
ax3.set_xlabel(r"x [$r_{\rm{g}}$]")
ax3.set_ylabel(r"y [$r_{\rm{g}}$]")
ax3.set_title(f"BLR emitting at rest frame {emitted_wavelength} nm \n observed in range {wavelength[0]}-{wavelength[1]} nm at z={redshift}")
ax3.set_aspect(1)

st.write(fig3)

st.write("Note that this figure may appear blank if the emission line does not project into your range of wavelengths.")

st.write("Here are the relevant R-Z projections of the density and poloidal velocities")

fig4, ax4 = plt.subplots(1, 2, figsize=(10, 4))

densities = my_blr.density_grid
velocities = np.sqrt(my_blr.z_velocity_grid**2 + my_blr.r_velocity_grid**2)

density_contours = ax4[0].contourf(R, Z, np.log10(densities*efficiency/np.nanmax(densities*efficiency)), 50)
velocity_contours = ax4[1].contourf(R, Z, velocities, 50)

density_cbar = plt.colorbar(density_contours, ax=ax4[0], label="log relative particle density * efficiency [arb.]")
velocity_cbar = plt.colorbar(velocity_contours, ax=ax4[1], label="outflowing velocities [c]")

plt.subplots_adjust(wspace=0.4)

for axis in ax4:
    axis.set_xlabel(r"R [$r_{\rm{g}}$]")
    axis.set_ylabel(r"Z [$r_{\rm{g}}$]")
    axis.set_aspect(1)

st.write(fig4)
    