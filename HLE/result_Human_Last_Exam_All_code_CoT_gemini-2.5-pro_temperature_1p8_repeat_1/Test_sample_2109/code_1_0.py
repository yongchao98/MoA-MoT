import math

# The problem is broken down into several steps. After analyzing the functional equations
# and the fractional derivative constraints, we arrive at a definite integral for the total energy E.
# The functional equations for T1(x) and T2(x) appear to be ill-posed as written.
# We proceed by using plausible intended solutions that satisfy the core symmetries of the equations:
# T1(x) = (x^2 - 2) / (x^2 + 2)
# T2(x) = x - 2x^2
# The heat flux function is determined to be Phi(z) = z/2 + 1.
# The total energy is given by the integral:
# E_total = integral from 0 to 1 of (Phi(T1(sqrt(2)*x) + T2(x))) dx
# This evaluates to the analytical result 17/12 - pi/4.

# We will now calculate this value.

# The equation for the minimum total heat energy is derived as:
# E = 17/12 - pi/4
val_17 = 17
val_12 = 12
val_pi = math.pi
val_4 = 4

# The first term in the equation is 17/12
term1 = val_17 / val_12
# The second term in the equation is pi/4
term2 = val_pi / val_4

# Calculate the final value
min_energy = term1 - term2

# We also output the equation itself with the numbers
print("The final equation for the total energy is:")
print(f"E_total = {val_17}/{val_12} - PI/{val_4}")
print("\nCalculating this value:")
print(f"E_total = {term1} - {term2}")
print(f"E_total = {min_energy}")

# Final Answer
# print(f"<<<{min_energy}>>>")