import numpy as np

# Step 1: Identify the poles from the E(z) field equation.
# Based on the analysis, the finite poles are z=1 and z=2.
poles_E = [1, 2]
sum_poles_E = sum(poles_E)
num_poles_E = len(poles_E)

# Step 2: Identify the poles from the B(z) field equation.
# The poles correspond to the roots of the polynomial 4*z^4 - z^3 + z^2 + 1 = 0.
# The coefficients are a_4=4, a_3=-1, a_2=1, a_1=0, a_0=1.
coeffs_B = [4, -1, 1, 0, 1]
num_poles_B = 4 # Degree of the polynomial

# Step 3: Calculate the sum of the poles from B(z) using Vieta's formulas.
# Sum of roots = -a_{n-1} / a_n
sum_poles_B = -coeffs_B[1] / coeffs_B[0]

# Step 4: Calculate the total number of poles and the total sum of poles.
total_num_poles = num_poles_E + num_poles_B
total_sum_poles = sum_poles_E + sum_poles_B

# Step 5: Compute the average value of the coordinates (poles).
average_value = total_sum_poles / total_num_poles

# Print the breakdown of the calculation as requested
print("The average value of the coordinates is calculated as follows:")
print(f"Sum of poles from E(z): {poles_E[0]} + {poles_E[1]} = {sum_poles_E}")
print(f"Sum of poles from B(z) (from Vieta's formulas): -({coeffs_B[1]}) / {coeffs_B[0]} = {sum_poles_B}")
print(f"Total sum of poles: {sum_poles_E} + {sum_poles_B} = {total_sum_poles}")
print(f"Total number of poles: {num_poles_E} + {num_poles_B} = {total_num_poles}")
print(f"Average value = (Total Sum) / (Total Number) = {total_sum_poles} / {total_num_poles} = {average_value}")

# Final Answer
# To provide a more exact fractional answer:
from fractions import Fraction
final_answer = Fraction(total_sum_poles).limit_denominator() / total_num_poles
print(f"The exact fractional value is: {final_answer}")
<<<13/24>>>