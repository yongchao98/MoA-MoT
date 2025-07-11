from fractions import Fraction

# --- Given values as fractions ---
# Density (rho) = 0.9 kg/cm^3
density = Fraction(9, 10)

# Radius (r) = 0.5 cm
radius = Fraction(1, 2)

# --- Constants and Approximations for the parrot ---
# The formula for the volume of a sphere uses 4/3
four_thirds = Fraction(4, 3)

# We need to approximate pi using integers <= 10.
# pi is approx 3.14159.
# Let's use pi_approx = 3.
# The error is |3 - pi| / pi, which is ~4.5%, less than the 10% limit.
pi_approx = Fraction(3, 1)

# --- Perform the calculation ---
# Volume = (4/3) * pi * r^3
# Mass = Density * Volume
# Let's combine them: Mass = density * (4/3) * pi_approx * radius^3
# The parrot calculates r^3 as (1/2) * (1/2) * (1/2)

# Calculate the final mass
mass = density * four_thirds * pi_approx * (radius**3)

# --- Output the equation as the parrot would calculate it ---
# The problem requires showing each number in the final equation.
print("The parrot's calculation for the mass is:")
print(f"{density.numerator}/{density.denominator} * {four_thirds.numerator}/{four_thirds.denominator} * {pi_approx.numerator} * ({radius.numerator}/{radius.denominator})^3 = {mass}")

# The integers involved are 9, 10, 4, 3, 3, 1, 2. The largest is 10.