import sympy

# Define all symbolic variables from the problem statement
K0, a, y, d, mu, mu0 = sympy.symbols('K_0 a y d mu mu_0', real=True, positive=True)

# The final formula for the force per unit area is constructed from several parts.
# 1. The constant pre-factor
pre_factor = mu0 / 2

# 2. The numerator of the main fraction
numerator = K0**2 * sympy.sin(a*y)**2

# 3. The denominator, which comes from applying the boundary conditions
denominator = (sympy.cosh(a*d) + (mu0/mu)*sympy.sinh(a*d))**2

# 4. The direction of the force is along the x-axis
direction_vector_str = "î_x"

# Assemble the final scalar part of the expression
force_per_area_scalar = pre_factor * (numerator / denominator)

# Print the final result in a clear, structured way
print("The final derived formula for the force per unit y-z area is:")
print("\n(Force / Area) =")
# Use sympy's pretty printer for a clean mathematical layout
sympy.pprint(force_per_area_scalar, use_unicode=True)
print(f"\nmultiplied by the direction vector: {direction_vector_str}")