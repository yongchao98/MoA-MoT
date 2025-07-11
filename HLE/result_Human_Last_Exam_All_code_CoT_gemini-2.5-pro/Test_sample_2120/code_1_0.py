import numpy as np

# Step 1: Identify the singularities from the electric field E(z).
# The singularities of E(z) are at z=1 and z=2.
singularities_E = [1, 2]
sum_E = sum(singularities_E)
num_E = len(singularities_E)

# Step 2: Identify the singularities from the magnetic field B(z).
# The singularities of B(z) are the roots of two polynomials.
# We don't need the roots, just their sum, which can be found using Vieta's formulas.
# Polynomial 1: 4*z^4 - 1*z^3 + 1*z^2 + 0*z + 1 = 0
# The sum of roots is -a_{n-1}/a_n = -(-1)/4 = 1/4.
poly1_coeffs = [4, -1, 1, 0, 1]
sum_roots_B1 = -poly1_coeffs[1] / poly1_coeffs[0]
num_roots_B1 = 4

# Polynomial 2: 1*z^4 + 0*z^3 + 1*z^2 - 1*z + 4 = 0
# The sum of roots is -a_{n-1}/a_n = -(0)/1 = 0.
poly2_coeffs = [1, 0, 1, -1, 4]
sum_roots_B2 = -poly2_coeffs[1] / poly2_coeffs[0]
num_roots_B2 = 4

# Step 3: Calculate the total sum and total number of singularities.
total_sum = sum_E + sum_roots_B1 + sum_roots_B2
total_num = num_E + num_roots_B1 + num_roots_B2

# Step 4: Compute the average value.
average_z = total_sum / total_num

# Step 5: Print the results, showing the numbers used in the final equation.
print("The singularities are located at:")
print(f"1. From E(z): {singularities_E[0]} and {singularities_E[1]}")
print(f"2. From B(z): the {num_roots_B1} roots of 4z^4 - z^3 + z^2 + 1 = 0")
print(f"3. From B(z): the {num_roots_B2} roots of z^4 + z^2 - z + 4 = 0")
print("\nCalculating the average value:")
print(f"Sum of singularities from E(z): {singularities_E[0]} + {singularities_E[1]} = {sum_E}")
print(f"Sum of the {num_roots_B1} singularities from the first polynomial for B(z): {sum_roots_B1}")
print(f"Sum of the {num_roots_B2} singularities from the second polynomial for B(z): {sum_roots_B2}")
print(f"Total number of singularities = {num_E} + {num_roots_B1} + {num_roots_B2} = {total_num}")
print(f"The final equation for the average value is: ({singularities_E[0]} + {singularities_E[1]} + {sum_roots_B1} + {sum_roots_B2}) / {total_num}")
print(f"Average value = {total_sum} / {total_num} = {average_z}")

# Final Answer
print("\n<<<" + str(average_z) + ">>>")