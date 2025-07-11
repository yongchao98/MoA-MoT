import sympy as sp

# Define symbolic variables for the pre-reset covariance matrix Sigma.
# For a covariance matrix, s_ij = s_ji.
s11, s12, s13 = sp.symbols('s11 s12 s13')
s22, s23 = sp.symbols('s22 s23')
s33 = sp.symbols('s33')

Sigma = sp.Matrix([
    [s11, s12, s13],
    [s12, s22, s23],
    [s13, s23, s33]
])

# As derived in the explanation, the estimation error itself is unchanged by the reset.
# This means the transformation on the covariance is an identity transformation.
# The post-reset covariance Sigma' is equal to the pre-reset covariance Sigma.
Sigma_prime = Sigma

# --- Output the result ---

print("The exact expression for the post-reset covariance Σ' is simply:")
print("Σ' = Σ")
print("\nThis means the covariance matrix is unchanged by the reset operation.")
print("To show the components of the equation, we represent the matrices symbolically:\n")

# Use symbolic matrix names for the final equation printout
Sigma_prime_symbol = sp.MatrixSymbol("Sigma_prime", 3, 3)
Sigma_symbol = sp.MatrixSymbol("Sigma", 3, 3)

# Print the final equation with matrix components
final_equation = sp.Eq(Sigma_prime_symbol, Sigma, evaluate=False)

print("Pre-reset covariance matrix:")
sp.pprint(sp.Eq(Sigma_symbol, Sigma))

print("\nPost-reset covariance matrix:")
sp.pprint(sp.Eq(Sigma_prime_symbol, Sigma_prime))
