import math

# Define the symbols used in the equation
pi = "pi"
d = "d"
h = "h"
rho = "rho"
t = "t"
g = "g"
H = "H"

# The numerical constants in the formula
four = 4
two = 2

# Construct and print the formula for the change in weight (Delta W)
# This formula corresponds to Option A.
# It represents the impact force from the sand landing.
# The mass flow rate is (pi * d**2 * h * rho) / (4 * t)
# The impact velocity is sqrt(2 * g * (H - h))

# Printing the numerator of the mass flow rate term
numerator = f"{pi} * {d}^2 * {h} * {rho}"

# Printing the denominator of the mass flow rate term
denominator = f"{four} * {t}"

# Printing the term inside the square root for velocity
sqrt_term = f"{two} * {g} * ({H} - {h})"

# Assembling and printing the final equation
print("The estimated change in weight is given by the expression:")
print(f"ΔW = ({numerator} / {denominator}) * sqrt({sqrt_term})")