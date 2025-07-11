import sympy

# Define 'n' as a symbolic variable to represent the dimension of the projective space
n = sympy.Symbol('n')

# The dimension of the space of global sections of O(d) on P^n is C(n+d, d).
def h0(d, var):
    """Calculates the dimension of H^0(P^n, O(d))."""
    return sympy.binomial(var + d, d)

# Dimension of the domain H^0(P^n, O(1)^{n+1})
h0_O1_n_plus_1 = (n + 1) * h0(1, n)

# Dimension of the codomain H^0(P^n, O(2))
h0_O2 = h0(2, n)

# The dimension of the space of global sections is the difference of the dimensions
# as derived from the long exact sequence in cohomology.
dim = h0_O1_n_plus_1 - h0_O2

# Print the calculation steps with symbolic formulas
print(f"The complex dimension is given by the formula h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^{{n+1}}) - h^0(P^n, O(2))")
print(f"First, calculate h^0(P^n, O(1)) = C(n+1, 1) = {h0(1, n)}")
print(f"Then, h^0(P^n, O(1)^{{n+1}}) = (n+1) * {h0(1, n)} = {sympy.simplify(h0_O1_n_plus_1)}")
print(f"Next, calculate h^0(P^n, O(2)) = C(n+2, 2) = {sympy.simplify(h0_O2)}")
print("The final dimension is the difference:")
print(f"dim = {sympy.simplify(h0_O1_n_plus_1)} - ({sympy.simplify(h0_O2)})")
print(f"dim = {sympy.simplify(dim)}")
print("\nThe final simplified formula for the dimension is:")
print(f"dim = n*(n+1)/2, or C(n+1, 2)")