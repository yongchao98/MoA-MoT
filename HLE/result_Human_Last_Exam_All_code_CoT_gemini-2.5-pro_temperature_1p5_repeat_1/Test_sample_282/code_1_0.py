import sympy

# Define symbolic variables for the dimensions of the set X
# Epsilon is a small positive value, L is a large positive value.
epsilon, L = sympy.symbols('epsilon L', positive=True, real=True)

# Our set X is a rectangle in the (y, x) coordinates for the AN subgroup of SL2(R).
# y is in [1-epsilon, 1+epsilon]
# x is in [0, L]
y_interval_X = sympy.Interval(1 - epsilon, 1 + epsilon)
x_interval_X = sympy.Interval(0, L)

# The measure of X is the area of this rectangle.
width_y_X = sympy.width(y_interval_X)
width_x_X = sympy.width(x_interval_X)
measure_X = width_y_X * width_x_X

print("Analysis of the set X:")
print(f"The y-interval of X is {y_interval_X}, with width {width_y_X}")
print(f"The x-interval of X is {x_interval_X}, with width {width_x_X}")
print(f"The measure of X is (width_y) * (width_x) = {width_y_X} * {width_x_X} = {measure_X}")
print("-" * 20)

# Now, we analyze the 3-fold product set X^3.
# An element (y', x') in X^3 is the product of three elements from X.
# y' = y1*y2*y3
# x' = y1*y2*x3 + y1/y3*x2 + 1/(y2*y3)*x1
# where yi are in y_interval_X and xi are in x_interval_X.

# For small epsilon, y values are close to 1.
# The y-interval for X^3 is approximately [(1-epsilon)^3, (1+epsilon)^3]
y_interval_X3 = sympy.Interval(y_interval_X.inf**3, y_interval_X.sup**3)
width_y_X3 = sympy.width(y_interval_X3)
# Let's expand and approximate for small epsilon
width_y_X3_approx = sympy.series(width_y_X3, epsilon, 0, 2).removeO()

# For the x-interval, we approximate y_i as 1.
# x' approx x1 + x2 + x3.
# The range of x' is approximately [0, 3L].
x_interval_X3_approx = sympy.Interval(0, 3*L)
width_x_X3_approx = sympy.width(x_interval_X3_approx)

# The measure of X^3 is approximated by the area of this new rectangle.
measure_X3_approx = width_y_X3_approx * width_x_X3_approx

print("Analysis of the set X^3 (approximated for small epsilon):")
print(f"The y-interval of X^3 is approximately {y_interval_X3}")
print(f"The width of the y-interval is {width_y_X3}, which is approximately {width_y_X3_approx}")
print(f"The x-interval of X^3 is approximately {x_interval_X3_approx}, with width {width_x_X3_approx}")
print(f"The approximate measure of X^3 is (width_y3) * (width_x3) = {width_y_X3_approx} * {width_x_X3_approx} = {measure_X3_approx}")
print("-" * 20)

# Finally, compute the ratio K.
K_approx = measure_X3_approx / measure_X
K_exact_limit = sympy.limit(K_approx, epsilon, 0)

print("Calculation of the constant K:")
print(f"The ratio mu(X^3)/mu(X) is approximately ({measure_X3_approx}) / ({measure_X})")
print(f"Simplifying the ratio gives: {sympy.simplify(K_approx)}")
print(f"The limit of this ratio as epsilon -> 0 is {K_exact_limit}")
print(f"\nFinal calculation steps:")
print(f"Measure of X = {measure_X}")
print(f"Approximate measure of X^3 = {measure_X3_approx}")
print(f"Ratio = ({sympy.expand(width_y_X3)}) * {width_x_X3_approx} / ({width_y_X} * {width_x_X})")
y_width_factor = sympy.limit(width_y_X3/width_y_X, epsilon, 0)
x_width_factor = sympy.limit(width_x_X3_approx/width_x_X, L, sympy.oo)
print(f"Ratio of y-interval widths in the limit: {y_width_factor}")
print(f"Ratio of x-interval widths in the limit: {x_width_factor}")
print(f"Total ratio K = {y_width_factor} * {x_width_factor}")
