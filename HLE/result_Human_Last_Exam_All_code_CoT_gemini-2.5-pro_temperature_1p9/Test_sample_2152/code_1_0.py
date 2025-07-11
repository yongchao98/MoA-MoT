import math

# Step 1: Define the given constants.
# The scattering angle in radians.
theta = math.pi / 30

# Step 2: Evaluate the dimensionless term (e_m * h_bar / (p * mu))^2.
# Let e be the elementary charge, h_bar be the reduced Planck constant,
# m be the particle's mass, and c be the speed of light.
# e_m = e / 16
# v = c / 100
# mu = 25 * mu_B = 25 * (e * h_bar / (2 * m * c)) (assuming particle is an electron)
# p = m * v = m * (c / 100)
#
# The term simplifies as:
# term = ( (e/16) * h_bar ) / ( (m*c/100) * (25 * e * h_bar / (2*m*c)) )
#      = (1/16) / ( (1/100) * (25/2) )
#      = (1/16) / (25/200)
#      = (1/16) / (1/8) = 8/16 = 0.5
#
# The term squared is:
dimless_term_sq = 0.5 ** 2

# Step 3: The ratio of differential cross-sections is given by the formula:
# Ratio = (3 / sin(theta)^2) * (dimensionless_term)^2
numerator = 3 * dimless_term_sq
denominator = math.sin(theta)**2
ratio = numerator / denominator

# Step 4: Print the final equation with each numerical component and the result.
print("The final equation is of the form: Numerator / Denominator = Ratio")
print(f"Where the numerator is 3 * (e_m*ħ / (p*μ))^2 and the denominator is sin(θ)^2.")
print("\nCalculating each part:")
print(f"Numerator = 3 * {dimless_term_sq} = {numerator}")
print(f"Denominator = sin({theta:.5f})^2 = {denominator:.6f}")
print("\nFinal ratio:")
print(f"{numerator} / {denominator:.6f} = {ratio:.4f}")
