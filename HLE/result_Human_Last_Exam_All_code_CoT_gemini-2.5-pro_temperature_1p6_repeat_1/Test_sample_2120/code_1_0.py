import numpy as np

# Define the coefficients of the polynomials P(z) and Q(z)
# P(z) = 4*z**4 - 1*z**3 + 1*z**2 + 0*z + 1
P_coeffs = [4, -1, 1, 0, 1]
# Q(z) = 1*z**4 + 0*z**3 + 1*z**2 - 1*z + 4
Q_coeffs = [1, 0, 1, -1, 4]

# Number of roots for each polynomial
num_roots_P = len(P_coeffs) - 1
num_roots_Q = len(Q_coeffs) - 1

# Calculate the sum of roots using Vieta's formulas (sum = -a_{n-1}/a_n)
# For P(z), a_n=4, a_{n-1}=-1
sum_roots_P = -P_coeffs[1] / P_coeffs[0]

# For Q(z), a_n=1, a_{n-1}=0
sum_roots_Q = -Q_coeffs[1] / Q_coeffs[0]

# The poles from E(z) are at z=1 and z=2
poles_E = [1, 2]
sum_poles_E = sum(poles_E)
num_poles_E = len(poles_E)

# Calculate the total sum and total number of poles
total_sum = sum_roots_P + sum_roots_Q + sum_poles_E
total_num_poles = num_roots_P + num_roots_Q + num_poles_E

# Calculate the average value
average_value = total_sum / total_num_poles

# Output the results of each step
print(f"Poles from B(z) analysis:")
print(f"Sum of the {num_roots_P} roots of P(z)=4z^4-z^3+z^2+1=0 is: {sum_roots_P}")
print(f"Sum of the {num_roots_Q} roots of Q(z)=z^4+z^2-z+4=0 is: {sum_roots_Q}")
print(f"Poles from E(z) analysis are at: {poles_E[0]} and {poles_E[1]}")
print(f"Total sum of poles = {sum_roots_P} + {sum_roots_Q} + {sum_poles_E} = {total_sum}")
print(f"Total number of poles = {num_roots_P} + {num_roots_Q} + {num_poles_E} = {total_num_poles}")
print(f"The average value of the coordinates z is {total_sum} / {total_num_poles} = {average_value}")

print(f"\nFinal Answer: {average_value}")
<<<0.325>>>