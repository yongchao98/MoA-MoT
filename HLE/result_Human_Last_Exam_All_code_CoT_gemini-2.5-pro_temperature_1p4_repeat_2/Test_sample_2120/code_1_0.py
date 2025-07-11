import numpy as np

# The poles of the B(z) field are determined by the roots of two polynomials,
# P1(z) and P2(z), derived from the denominators of the functions R(z) and R(1/z).

# P1(z) = 4*z**4 - z**3 + z**2 + 1
P1_coeffs = [4, -1, 1, 0, 1]

# P2(z) = z**4 + z**2 - z + 4
P2_coeffs = [1, 0, 1, -1, 4]

# According to Vieta's formulas, the sum of the roots of a polynomial
# a_n*z**n + a_{n-1}*z**(n-1) + ... + a_0 is -a_{n-1} / a_n.

# Calculate the sum of roots for P1(z)
sum_roots_P1 = -P1_coeffs[1] / P1_coeffs[0]

# Calculate the sum of roots for P2(z)
sum_roots_P2 = -P2_coeffs[1] / P2_coeffs[0]

# The total number of poles is the sum of the degrees of the two polynomials.
num_poles_P1 = len(P1_coeffs) - 1
num_poles_P2 = len(P2_coeffs) - 1
total_num_poles = num_poles_P1 + num_poles_P2

# The total sum of poles is the sum of the sums of roots.
total_sum_poles = sum_roots_P1 + sum_roots_P2

# The average value is the total sum divided by the total number of poles.
average_value = total_sum_poles / total_num_poles

# Output the step-by-step calculation
print("The poles are the roots of two polynomials, P1(z) and P2(z).")
print("P1(z) = 4*z^4 - 1*z^3 + 1*z^2 + 0*z + 1")
print("P2(z) = 1*z^4 + 0*z^3 + 1*z^2 - 1*z + 4")
print("\nUsing Vieta's formulas to find the sum of roots for each polynomial:")

print("\nSum of roots for P1(z) = - (coefficient of z^3) / (coefficient of z^4)")
print(f"Sum of roots for P1(z) = -({P1_coeffs[1]}) / ({P1_coeffs[0]}) = {sum_roots_P1}")

print("\nSum of roots for P2(z) = - (coefficient of z^3) / (coefficient of z^4)")
print(f"Sum of roots for P2(z) = -({P2_coeffs[1]}) / ({P2_coeffs[0]}) = {sum_roots_P2}")

print("\nTotal sum of poles = (Sum of roots of P1) + (Sum of roots of P2)")
print(f"Total sum of poles = {sum_roots_P1} + {sum_roots_P2} = {total_sum_poles}")

print("\nTotal number of poles = (Degree of P1) + (Degree of P2)")
print(f"Total number of poles = {num_poles_P1} + {num_poles_P2} = {total_num_poles}")

print("\nAverage value = Total sum of poles / Total number of poles")
print(f"Average value = ({total_sum_poles}) / ({total_num_poles})")
print(f"Average value = {average_value}")