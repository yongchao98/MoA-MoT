import sympy as sp

# Step 1: Define the symbolic variables.
# u: the coupling constant of the φ⁴ interaction.
# c1: a non-zero constant representing the one-loop contribution coefficient.
u = sp.Symbol('u')
c1 = sp.Symbol('c1')

# Step 2: Define the expression for the critical exponent ν as a function of u.
# The exponent ν is the inverse of the thermal scaling dimension y_t = 2 + η_φ²(u).
# The anomalous dimension η_φ²(u) starts at linear order in u, so we write η_φ²(u) ≈ c1*u.
eta_phi2 = c1 * u
nu_expr = 1 / (2 + eta_phi2)

# Step 3: Perform a Taylor series expansion for ν(u) around u=0.
# This reveals the structure of the corrections to the mean-field value.
nu_series = nu_expr.series(u, 0, 2)

# Step 4: Extract the terms from the series to construct and print the final equation.
zeroth_order_term = nu_series.coeff(u, 0)
first_order_coeff = nu_series.coeff(u, 1)

# Step 5: Print the result and explanation.
print("The critical exponent ν as a function of the coupling constant u is given by ν(u) = 1 / (2 + η_φ²(u)).")
print("The anomalous dimension η_φ²(u) is a power series in u, starting at order u¹:")
print(f"η_φ²(u) = {eta_phi2} + O(u²)")
print("\nExpanding ν(u) as a power series for small u results in the following equation:")
print("\nFinal Equation:")
print(f"ν = {zeroth_order_term} + ({first_order_coeff})*u + O(u²)")
print("\nThis equation shows that the mean-field value is ν = 1/2.")
print(f"The initial non-vanishing contribution beyond this value is the term ({first_order_coeff})*u.")
print("The order of this term with respect to the coupling constant u is 1.")
