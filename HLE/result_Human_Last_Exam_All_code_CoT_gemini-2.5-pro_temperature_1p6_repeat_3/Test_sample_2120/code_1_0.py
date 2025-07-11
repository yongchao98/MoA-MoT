import numpy as np

# Step 1: Define the problem as finding the average of the poles of S(z) = E(z)B(z)

# Step 2: Analyze B(z) to find the sum and number of its poles.
# The poles of B(z) are the roots of two polynomials.
# P1(z) = 4*z^4 - z^3 + z^2 + 1
poly_B1_coeffs = [4, -1, 1, 0, 1]
# P2(z) = z^4 + z^2 - z + 4
poly_B2_coeffs = [1, 0, 1, -1, 4]

# According to Vieta's formulas, the sum of the roots of a polynomial
# a_n*z^n + a_{n-1}*z^{n-1} + ... + a_0 is -a_{n-1}/a_n.
# The degree of the polynomial is n.

# For P1(z)
deg_B1 = len(poly_B1_coeffs) - 1
sum_poles_B1 = -poly_B1_coeffs[1] / poly_B1_coeffs[0]
num_poles_B1 = deg_B1

# For P2(z)
deg_B2 = len(poly_B2_coeffs) - 1
sum_poles_B2 = -poly_B2_coeffs[1] / poly_B2_coeffs[0]
num_poles_B2 = deg_B2

# Total sum and number of poles for B(z)
sum_poles_B = sum_poles_B1 + sum_poles_B2
num_poles_B = num_poles_B1 + num_poles_B2

# Step 3: Analyze E(z) to find the sum and number of its poles.
# From the analysis of the functional equation for E(z), its poles are at z=1 and z=2.
poles_E = [1, 2]
sum_poles_E = sum(poles_E)
num_poles_E = len(poles_E)

# Step 4: Calculate the total average.
# The sets of poles are disjoint.
total_sum_of_poles = sum_poles_E + sum_poles_B
total_num_of_poles = num_poles_E + num_poles_B

# Calculate the average value
average_value = total_sum_of_poles / total_num_of_poles

# Output the results, showing each number in the final equation.
print("Calculating the average value of the complex coordinates z.")
print("This is the sum of all poles divided by the number of poles.")
print("\nPoles from E(z):")
print(f"  Sum: {poles_E[0]} + {poles_E[1]} = {sum_poles_E}")
print(f"  Number of poles: {num_poles_E}")

print("\nPoles from B(z):")
# Showing the fractions for clarity
sum_B1_str = f"({-poly_B1_coeffs[1]}/{poly_B1_coeffs[0]})"
sum_B2_str = f"({-poly_B2_coeffs[1]}/{poly_B2_coeffs[0]})"
print(f"  Sum: {sum_B1_str} + {sum_B2_str} = {sum_poles_B1} + {sum_poles_B2} = {sum_poles_B}")
print(f"  Number of poles: {num_poles_B1} + {num_poles_B2} = {num_poles_B}")

print("\nTotal calculation:")
print(f"  Total Sum = {sum_poles_E} + {sum_poles_B} = {total_sum_of_poles}")
print(f"  Total Number = {num_poles_E} + {num_poles_B} = {total_num_of_poles}")
print(f"  Average = ({sum_poles_E} + {sum_poles_B1} + {sum_poles_B2}) / ({num_poles_E} + {num_poles_B}) = {total_sum_of_poles} / {total_num_of_poles}")

print(f"\nThe final average value is: {average_value}")

print("\nFinal Equation:")
print(f"({sum_poles_E} + {sum_poles_B}) / {total_num_of_poles} = {average_value}")