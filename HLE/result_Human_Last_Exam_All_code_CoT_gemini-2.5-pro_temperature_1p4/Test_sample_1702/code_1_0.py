import math

# Define the physical constants and given values
# h: Planck's constant in J·s
# c: Speed of light in m/s
# k: Boltzmann constant in J/K
# B: Spectral radiance in W/m^2/sr/m (which is W·m⁻³·sr⁻¹)
# lambda_: Wavelength in m
h = 6.62607015e-34
c = 299792458.0
k = 1.380649e-23
B = 9.9e16
lambda_ = 500e-9 # 500 nm

# To solve for Temperature (T), we use the rearranged Planck's Law:
# T = (h * c / (lambda_ * k)) / ln((2 * h * c**2) / (B * lambda_**5) + 1)
# Let's calculate the numerator and denominator of this final division separately.

# Calculate the numerator part: (h * c) / (lambda_ * k)
numerator = (h * c) / (lambda_ * k)

# Calculate the argument for the natural logarithm: (2 * h * c**2) / (B * lambda_**5) + 1
log_argument = (2 * h * c**2) / (B * lambda_**5) + 1

# Calculate the denominator part: ln(log_argument)
denominator = math.log(log_argument)

# Now, calculate the temperature in Kelvin by dividing the numerator by the denominator
temperature_in_kelvin = numerator / denominator

# Convert the result to thousands of Kelvin and round to the nearest whole number
final_answer = round(temperature_in_kelvin / 1000)

# As requested, output each number in the final equation for temperature calculation
print("The final equation to calculate temperature T is T = A / B")
print(f"Where A = (h * c / (λ * k)) = {numerator}")
print(f"And   B = ln((2 * h * c^2) / (B * λ^5) + 1) = {denominator}")
print(f"Resulting Temperature (Kelvin) = {temperature_in_kelvin}\n")

# Print the final answer in the requested format
print("The temperature of Pandora in thousands of Kelvin (rounded) is:")
print(final_answer)