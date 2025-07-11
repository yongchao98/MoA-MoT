from decimal import Decimal, getcontext, ROUND_HALF_UP

# Set precision for Decimal calculations to ensure accuracy.
getcontext().prec = 50

# Problem parameters and physical constants.
# The value for k is approximated to 1.4e-23 to align with the
# representational constraints of the fictional 'frac' data type.
B = Decimal("9.9e16")  # Spectral radiance in W/m^2/sr/m
L = Decimal("5e-7")    # Wavelength in meters (500 nm)
C = Decimal("3e8")     # Speed of light in m/s
K = Decimal("1.4e-23") # Boltzmann constant in J/K (approximated)
TWO = Decimal("2")
THOUSAND = Decimal("1000")

# Calculate temperature in Kelvin using the Rayleigh-Jeans Law:
# T = (B * L^4) / (2 * C * K)
numerator = B * (L ** 4)
denominator = TWO * C * K
temp_kelvin = numerator / denominator

# Convert temperature to thousands of Kelvin.
temp_kilo_kelvin = temp_kelvin / THOUSAND

# Round the result to the nearest integer.
rounded_temp_kk = temp_kilo_kelvin.to_integral_value(rounding=ROUND_HALF_UP)

# Per the instructions, output the numbers in the final equation along with the result.
# We format the constants as strings for clear presentation in the equation.
B_str = "9.9e16"
L_str = "5e-7"
C_str = "3e8"
K_str = "1.4e-23"
TWO_str = "2"
THOUSAND_str = "1000"

print(f"round( (({B_str}) * ({L_str})^4) / (({TWO_str}) * ({C_str}) * ({K_str})) / {THOUSAND_str} ) = {rounded_temp_kk}")