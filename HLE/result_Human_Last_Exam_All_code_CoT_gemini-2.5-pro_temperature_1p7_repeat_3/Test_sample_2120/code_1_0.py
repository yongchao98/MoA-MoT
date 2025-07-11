import numpy as np

# --- Plan ---
# 1. Identify the singularity coordinates from the equations for E(z) and B(z).
# 2. For E(z), the singularities are z=1 and z=2.
# 3. For B(z), the singularities are the roots of two polynomials derived from the RHS of the equation.
# 4. Use Vieta's formulas to find the sum of the roots without solving for them.
# 5. Sum all singularity coordinates and divide by the total number of singularities.

# --- Calculations ---

# Singularities from E(z)
s_e = [1, 2]
sum_e = sum(s_e)
count_e = len(s_e)
print(f"Singularities from E(z) are {s_e[0]} and {s_e[1]}.")
print(f"Sum of singularities from E(z): {s_e[0]} + {s_e[1]} = {sum_e}")
print(f"Number of singularities from E(z): {count_e}")
print("-" * 20)

# Singularities from B(z)
# These are the roots of P1(z) = 4*z^4 - z^3 + z^2 + 1 = 0
# and their reciprocals, which are the roots of P2(z) = z^4 + z^2 - z + 4 = 0.

# For P1(z) = 4*z^4 - z^3 + z^2 + 0*z + 1 = 0
poly1_coeffs = [4, -1, 1, 0, 1]
# Sum of roots = -a_{n-1}/a_n
sum_b1 = -poly1_coeffs[1] / poly1_coeffs[0]
count_b1 = 4
print(f"The first set of B(z) singularities are the roots of 4z^4 - z^3 + z^2 + 1 = 0.")
print(f"Sum of these roots (by Vieta's formulas): -(-1) / 4 = {sum_b1}")
print(f"Number of these roots: {count_b1}")
print("-" * 20)

# For P2(z) = z^4 + 0*z^3 + z^2 - z + 4 = 0
poly2_coeffs = [1, 0, 1, -1, 4]
# Sum of roots = -a_{n-1}/a_n
sum_b2 = -poly2_coeffs[1] / poly2_coeffs[0]
count_b2 = 4
print(f"The second set of B(z) singularities are the roots of z^4 + z^2 - z + 4 = 0.")
print(f"Sum of these roots (by Vieta's formulas): -(0) / 1 = {sum_b2}")
print(f"Number of these roots: {count_b2}")
print("-" * 20)

# --- Final Calculation ---

# Total sum of singularities
total_sum = sum_e + sum_b1 + sum_b2
print("Calculating total sum of coordinates:")
print(f"Total Sum = (Sum from E) + (Sum from B set 1) + (Sum from B set 2)")
print(f"Total Sum = {sum_e} + {sum_b1} + {sum_b2} = {total_sum}")
print("-" * 20)

# Total number of singularities
total_count = count_e + count_b1 + count_b2
print("Calculating total number of coordinates:")
print(f"Total Count = (Count from E) + (Count from B set 1) + (Count from B set 2)")
print(f"Total Count = {count_e} + {count_b1} + {count_b2} = {total_count}")
print("-" * 20)

# Average value
average = total_sum / total_count
print("Calculating the average value:")
print(f"Average = Total Sum / Total Count")
print(f"Average = {total_sum} / {total_count} = {average}")

print("\nFinal Answer:")
print(f"The final equation for the average is: (({s_e[0]} + {s_e[1]}) + ({sum_b1}) + ({sum_b2})) / {total_count}")
print(f"Which simplifies to: {total_sum} / {total_count}")
<<<0.325>>>