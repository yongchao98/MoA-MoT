import numpy as np

# According to the plan, we calculate the average of the coordinates of the
# singularities of B(z), which are the roots of two polynomials.
# P1(z) = 4z^4 - z^3 + z^2 + 1
# P2(z) = z^4 + z^2 - z + 4

# We use Vieta's formulas to find the sum of roots for each polynomial.
# For a polynomial a_n*z^n + a_{n-1}*z^{n-1} + ... + a_0, the sum of roots is -a_{n-1}/a_n.

# For P1(z): a_4=4, a_3=-1. The polynomial has 4 roots.
sum_roots_P1 = -(-1) / 4
num_roots_P1 = 4

# For P2(z): a_4=1, a_3=0. The polynomial has 4 roots.
sum_roots_P2 = -(0) / 1
num_roots_P2 = 4

# The total sum of the coordinates of all singularities.
total_sum = sum_roots_P1 + sum_roots_P2

# The total number of singularities.
total_num_roots = num_roots_P1 + num_roots_P2

# The average value of the coordinates.
average_value = total_sum / total_num_roots

print("The final calculation for the average value is based on the sums of roots of the derived polynomials.")
print(f"Sum of roots from P1(z): {sum_roots_P1}")
print(f"Number of roots from P1(z): {num_roots_P1}")
print(f"Sum of roots from P2(z): {sum_roots_P2}")
print(f"Number of roots from P2(z): {num_roots_P2}")
print("Final equation for the average:")
print(f"({sum_roots_P1} + {sum_roots_P2}) / ({num_roots_P1} + {num_roots_P2}) = {average_value}")

print("<<<0.03125>>>")