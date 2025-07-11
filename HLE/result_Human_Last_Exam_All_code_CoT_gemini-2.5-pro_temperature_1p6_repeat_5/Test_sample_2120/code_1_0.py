# Plan:
# 1. Identify the number and sum of poles from the E(z) field equation.
# 2. Identify the number and sum of poles from the B(z) field equation using Vieta's formulas.
# 3. Calculate the total number and total sum of all poles.
# 4. Compute the average value.
# 5. Print the final equation with all its components.

# Poles from E(z)
poles_E = [1, 2]
num_poles_E = len(poles_E)
sum_poles_E = sum(poles_E)

# Poles from B(z)
# The poles are the roots of two polynomials.
# Polynomial 1: 4*z**4 - z**3 + z**2 + 1 = 0
# From Vieta's formulas, for a polynomial a*z**n + b*z**(n-1) + ...,
# the sum of roots is -b/a.
a1, b1 = 4, -1
sum_roots_P1 = -b1 / a1
num_roots_P1 = 4

# Polynomial 2: z**4 + z**2 - z + 4 = 0
a2, b2 = 1, 0
sum_roots_P2 = -b2 / a2
num_roots_P2 = 4

num_poles_B = num_roots_P1 + num_roots_P2
sum_poles_B = sum_roots_P1 + sum_roots_P2

# Total poles
total_num_poles = num_poles_E + num_poles_B
total_sum_poles = sum_poles_E + sum_poles_B

# Calculate the average
average = total_sum_poles / total_num_poles

# Output the results clearly, showing the final equation
print("The average of the complex coordinates is calculated by summing the poles and dividing by the number of poles.")
print(f"Sum of poles from E(z): {poles_E[0]} + {poles_E[1]} = {sum_poles_E}")
print(f"Sum of poles from B(z): ({sum_roots_P1}) + ({sum_roots_P2}) = {sum_poles_B}")
print(f"Total sum of poles: {sum_poles_E} + {sum_poles_B} = {total_sum_poles}")
print(f"Total number of poles: {num_poles_E} + {num_poles_B} = {total_num_poles}")
print("-" * 20)
print("Final Calculation:")
# Display the full calculation as a single equation
print(f"({poles_E[0]} + {poles_E[1]} + {sum_roots_P1} + {sum_roots_P2}) / ({num_poles_E} + {num_roots_P1} + {num_roots_P2}) = {total_sum_poles} / {total_num_poles} = {average}")

<<<0.325>>>