import sympy

# Define the variables for our polynomial over F_2
X, Y = sympy.symbols('X Y')

# Define the polynomial for the reduced curve
# In F_2, we work with integers mod 2. 5x^4 becomes x^4, etc.
F = Y**2 - (X**5 + X**4 + X**2)

# Calculate the partial derivatives
F_X = sympy.diff(F, X)
F_Y = sympy.diff(F, Y)

# In F_2, we perform calculations modulo 2
# sympy's Poly can handle this
poly_F = sympy.Poly(F, X, Y, modulus=2)
poly_F_X = sympy.Poly(F_X, X, Y, modulus=2)
poly_F_Y = sympy.Poly(F_Y, X, Y, modulus=2)

# Find the singular points by solving the system {F=0, dF/dX=0, dF/dY=0} mod 2
singular_points = sympy.solve([poly_F, poly_F_X, poly_F_Y], [X, Y])

print(f"The equation of the reduced curve is: Y^2 = X^5 + X^4 + X^2")
print(f"Partial derivative w.r.t X (mod 2): {poly_F_X.as_expr()}")
print(f"Partial derivative w.r.t Y (mod 2): {poly_F_Y.as_expr()}")
print(f"The affine singular points are: {singular_points}")

# From our step-by-step analysis, we know:
# 1. There is one affine singularity at (0,0) which resolves to one double point.
# 2. There is one singularity at infinity which also resolves to one double point.
num_affine_nodes = len(singular_points)
num_infinite_nodes = 1 # From manual analysis
total_nodes = num_affine_nodes + num_infinite_nodes

print(f"\nThe affine singularity at (0,0) contributes 1 double point.")
print(f"The singularity at infinity contributes 1 double point.")
print(f"Total number of double points = {num_affine_nodes} + {num_infinite_nodes} = {total_nodes}")
