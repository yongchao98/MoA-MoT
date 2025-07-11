import math

# 1. Define physical constants and problem givens
h = 6.626e-34  # Planck's constant in J·s
c = 2.998e8   # Speed of light in m/s
k = 1.381e-23  # Boltzmann constant in J/K
lambda_val = 500e-9 # Wavelength in meters (500 nm)
B_obs = 9.9e16      # Observed spectral radiance in W/m^2·sr·m

# 2. Calculate the intermediate components based on the rearranged Planck's Law
# Let C1 = (2*h*c^2 / λ^5)
C1 = (2 * h * c**2) / (lambda_val**5)
# Let C2 = (h*c / λ*k)
C2 = (h * c) / (lambda_val * k)

# 3. Calculate the argument for the natural logarithm function
# Z = C1 / B_obs
Z = C1 / B_obs

# 4. Calculate ln(1 + Z) using its Taylor Series expansion.
# Using the first 5 terms for good accuracy.
ln_val = Z - (Z**2)/2 + (Z**3)/3 - (Z**4)/4 + (Z**5)/5

# 5. Calculate the final temperature in Kelvin
T_kelvin = C2 / ln_val

# 6. Output the numbers used in the final step of the calculation
# The final equation is T = C2 / ln(1 + Z)
print("The final calculation is performed using the formula: T = C2 / ln(1 + Z)")
print("Final equation with numerical values:")
print(f"{T_kelvin} = {C2} / {ln_val}")
print("") # Adding a blank line for clarity

# 7. Convert the temperature to thousands of Kelvin and round to the nearest integer
final_answer = round(T_kelvin / 1000)

print("The final temperature of Pandora in thousands of Kelvin (rounded) is:")
print(final_answer)