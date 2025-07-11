import math

# The problem asks for the parameters (R, φ) of a deferent-epicycle model
# approximating an object moving on a square.
# R is the ratio of the deferent radius to the epicycle radius.
# φ is the ratio of the epicycle frequency to the deferent frequency.

# The parameter R can be determined from the ratio of the maximum and minimum
# distance of the object from the center of the square.
# Let the square have a side length of 2, centered at the origin.
# The maximum distance from the center (to a vertex) is sqrt(1^2 + 1^2) = sqrt(2).
# The minimum distance from the center (to a side midpoint) is 1.
# The ratio r_max / r_min = sqrt(2).

# In the epicycle model, r_max / r_min = (R + 1) / (R - 1), where R is the ratio of the radii.
# By setting these equal, we get the equation for R:
# (R + 1) / (R - 1) = sqrt(2)
# Solving for R:
# R + 1 = sqrt(2) * R - sqrt(2)
# R * (sqrt(2) - 1) = sqrt(2) + 1
# R = (sqrt(2) + 1) / (sqrt(2) - 1) = (sqrt(2) + 1)^2 = 3 + 2*sqrt(2).

sqrt_2 = math.sqrt(2)
R_val = 3 + 2 * sqrt_2

# The parameter φ is determined from the frequency of the radial oscillation.
# The object's distance from the center completes 4 cycles per full orbit.
# Thus, the radial frequency ω_r is 4 times the deferent frequency ω_d.
# In the model, the radial frequency is given by ω_r = |ω_e - ω_d| = |φ - 1| * ω_d.
# This gives the equation for φ: |φ - 1| = 4.
# This yields φ = 5 or φ = -3. The retrograde solution (φ < 0) is chosen as it's
# needed to approximate the sharp corners of the square.

phi_val = -3.0

# Print the final equations and the results.
print("1. Calculation of R:")
print(f"The equation for R derived from the geometry is: (R + 1) / (R - 1) = sqrt(2)")
print(f"Solving this gives R = (sqrt(2) + 1) / (sqrt(2) - 1), which simplifies to:")
print(f"R = 3 + 2 * {sqrt_2}")
print(f"The numerical value is R = {R_val}\n")

print("2. Calculation of φ:")
print(f"The equation for φ derived from the geometry is: |φ - 1| = 4")
print("Choosing the retrograde solution (φ < 0) gives:")
print(f"φ = {phi_val}\n")

print("3. The ordered pair (R, φ) is:")
print(f"({R_val}, {phi_val})")