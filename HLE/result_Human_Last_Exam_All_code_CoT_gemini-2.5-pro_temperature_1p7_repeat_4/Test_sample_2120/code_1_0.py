import numpy as np

# Step 1: Define the coefficients of the polynomial from the denominator.
# The polynomial is P(z) = 4*z^4 - z^3 + z^2 + 1 = 0.
a_4 = 4  # coefficient of z^4
a_3 = -1 # coefficient of z^3
n = 4    # degree of the polynomial

# Step 2: Calculate the sum of the roots using Vieta's formulas.
# Sum of roots = -(a_{n-1} / a_n)
sum_of_roots = -(a_3 / a_4)

# Step 3: Calculate the average of the roots.
# Average = (Sum of roots) / n
average_of_roots = sum_of_roots / n

# Step 4: Print the reasoning and the final calculation.
print("The coordinates 'z' are the poles of the fields, determined by the denominator of the source term in the B-field equation.")
print("This gives the polynomial equation: 4*z^4 - z^3 + z^2 + 1 = 0.")
print("According to Vieta's formulas, the sum of the roots (z1+z2+z3+z4) is -(-1 / 4).")
print(f"Sum of roots = {sum_of_roots}")
print("The average of the roots is the sum divided by the number of roots, which is 4.")
print("\nThe final equation for the average value is:")
print(f"{sum_of_roots} / {n} = {average_of_roots}")

# Output the final numerical answer in the specified format
# The output is formatted as required.