import numpy as np

# Step 1: Define the coefficients of the two polynomials determining the poles of B(z).
# P1(z) = 4*z^4 - 1*z^3 + 1*z^2 + 0*z + 1
p1_coeffs = np.array([4, -1, 1, 0, 1])

# P2(z) = 1*z^4 + 0*z^3 + 1*z^2 - 1*z + 4
p2_coeffs = np.array([1, 0, 1, -1, 4])

# Step 2: Use Vieta's formulas to find the sum of the roots for each polynomial.
# For a polynomial a_n*z^n + a_{n-1}*z^{n-1} + ..., the sum of roots is -a_{n-1}/a_n.

# Sum of roots for P1(z)
# a_4 = 4, a_3 = -1
sum_roots_p1_num = -p1_coeffs[1]
sum_roots_p1_den = p1_coeffs[0]
sum_roots_p1 = sum_roots_p1_num / sum_roots_p1_den

# Sum of roots for P2(z)
# a_4 = 1, a_3 = 0
sum_roots_p2_num = -p2_coeffs[1]
sum_roots_p2_den = p2_coeffs[0]
sum_roots_p2 = sum_roots_p2_num / sum_roots_p2_den

# Step 3: Calculate the total sum of the pole coordinates.
# We assume the poles of E(z) do not contribute additional singularities.
total_sum_of_roots = sum_roots_p1 + sum_roots_p2

# Step 4: Calculate the total number of poles.
num_poles = len(p1_coeffs) - 1 + len(p2_coeffs) - 1

# Step 5: Compute the average value of the pole coordinates.
average_value = total_sum_of_roots / num_poles

# Output the steps of the calculation
print("The poles of the B-field are determined by the roots of two polynomials:")
print(f"P1(z) = {p1_coeffs[0]}z^4 + {p1_coeffs[1]}z^3 + {p1_coeffs[2]}z^2 + {p1_coeffs[3]}z + {p1_coeffs[4]} = 0")
print(f"P2(z) = {p2_coeffs[0]}z^4 + {p2_coeffs[1]}z^3 + {p2_coeffs[2]}z^2 + {p2_coeffs[3]}z + {p2_coeffs[4]} = 0")
print("")
print("Using Vieta's formulas:")
print(f"Sum of roots for P1(z) = -({p1_coeffs[1]}) / {p1_coeffs[0]} = {sum_roots_p1}")
print(f"Sum of roots for P2(z) = -({p2_coeffs[1]}) / {p2_coeffs[0]} = {sum_roots_p2}")
print("")
print(f"Total sum of pole coordinates = {sum_roots_p1} + {sum_roots_p2} = {total_sum_of_roots}")
print(f"Total number of poles = 4 + 4 = {num_poles}")
print("")
print(f"The average value of the complex coordinates is the total sum divided by the number of poles:")
print(f"Average = {total_sum_of_roots} / {num_poles} = {average_value}")

print(f"\nThe final computed average value is {average_value}")