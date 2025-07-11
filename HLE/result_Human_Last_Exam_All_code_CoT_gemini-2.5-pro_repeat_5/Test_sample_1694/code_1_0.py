import cmath
import math

# --- Equation 1 ---
# (5 0; 0 -5)X1^2 + X1^2(6 0; 0 6) = (-53/12 0; 0 0)
# This simplifies to a linear equation for the elements of Y1 = X1^2.
# For the (1,1) element y1_11, we have: 5*y1_11 + 6*y1_11 = -53/12
# 11*y1_11 = -53/12
y1_11 = (-53/12) / 11

# The first coordinate of the solutions for X1 are the square roots of y1_11.
x1_sol1 = cmath.sqrt(y1_11)
x1_sol2 = -x1_sol1

print(f"For the first equation:")
print(f"The (1,1) element of X1^2 is {y1_11:.4f}")
print(f"The two possible first coordinates for X1 are {x1_sol1:.4f} and {x1_sol2:.4f}")
print("-" * 20)

# --- Equation 2 ---
# (4 0; 0 -5)X2^2 + X2^2(6 0; 0 6) = (-3/11 0; 0 0)
# For the (1,1) element y2_11, we have: 4*y2_11 + 6*y2_11 = -3/11
# 10*y2_11 = -3/11
y2_11 = (-3/11) / 10

# The first coordinate of the solutions for X2 are the square roots of y2_11.
x2_sol1 = cmath.sqrt(y2_11)
x2_sol2 = -x2_sol1

print(f"For the second equation:")
print(f"The (1,1) element of X2^2 is {y2_11:.4f}")
print(f"The two possible first coordinates for X2 are {x2_sol1:.4f} and {x2_sol2:.4f}")
print("-" * 20)

# --- Summation ---
# The total sum is the sum of all possible first coordinates found.
total_sum = x1_sol1 + x1_sol2 + x2_sol1 + x2_sol2

print("The final equation for the sum is:")
print(f"Sum = ({x1_sol1:.4f}) + ({x1_sol2:.4f}) + ({x2_sol1:.4f}) + ({x2_sol2:.4f})")
# We use .real and .imag to avoid potential floating point noise like -0.0j
print(f"Sum = {total_sum.real:.4f}{total_sum.imag:+.4f}j")

<<<0>>>