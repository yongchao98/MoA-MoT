import sympy

# Define n as a symbolic variable, assuming n is a positive integer.
n = sympy.Symbol('n', integer=True, positive=True)

# The dimension of the space of global sections is given by the formula:
# h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))

# The dimension of H^0(P^n, O(k)) is binomial(n + k, n) for k >= 0.

# 1. Calculate the dimension of the first term: h^0(P^n, O(1)^(n+1))
# h^0(P^n, O(1)) = binomial(n + 1, n)
h0_O1 = sympy.binomial(n + 1, n)
h0_O1_sum = (n + 1) * h0_O1
h0_O1_sum_simplified = sympy.simplify(h0_O1_sum)

# 2. Calculate the dimension of the second term: h^0(P^n, O(2))
h0_O2 = sympy.binomial(n + 2, n)
h0_O2_simplified = sympy.simplify(h0_O2)

# 3. The final dimension is the difference.
dim = h0_O1_sum - h0_O2
final_dim = sympy.simplify(dim)

# Print the components of the final equation and the result.
# The "numbers" in the equation are the dimensions of the two spaces.
print(f"The dimension of H^0(P^n, O(1)^(n+1)) is: {h0_O1_sum_simplified}")
print(f"The dimension of H^0(P^n, O(2)) is: {h0_O2_simplified}")
print(f"The final equation for the dimension is: {h0_O1_sum_simplified} - ({h0_O2_simplified})")
print(f"The complex dimension h^0(P^n, Omega^1(2)) simplifies to: {final_dim}")
