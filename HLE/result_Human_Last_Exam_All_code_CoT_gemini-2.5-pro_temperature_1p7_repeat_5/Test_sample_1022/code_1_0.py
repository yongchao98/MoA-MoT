# The degree of the bounded cohomology group we want to compute for T x T.
n = 4

# Known dimensions of the bounded cohomology of Thompson's group T with real coefficients.
# dim H_b^k(T; R) is stored in a dictionary. For k >= 4, the dimension is 0.
dim_H_b_T = {
    0: 1,
    1: 0,
    2: 1,
    3: 1,
    4: 0,
}

# We apply the Künneth formula for bounded cohomology:
# dim H_b^n(T x T; R) = Σ_{p+q=n} dim H_b^p(T; R) * dim H_b^q(T; R)

total_dim = 0
equation_parts = []
symbolic_parts = []

# Iterate through all non-negative integers p from 0 to n.
for p in range(n + 1):
    q = n - p
    
    # Retrieve the dimensions for degree p and q.
    # The .get(key, 0) method returns 0 if the key is not found.
    dim_p = dim_H_b_T.get(p, 0)
    dim_q = dim_H_b_T.get(q, 0)
    
    # Calculate the product for the current term and add to the total dimension.
    term_product = dim_p * dim_q
    total_dim += term_product
    
    # Build up the string representation of the calculation.
    equation_parts.append(f"({dim_p} * {dim_q})")
    symbolic_parts.append(f"dim H_b^{p}(T) * dim H_b^{q}(T)")

# Join the parts to form the full equation string.
full_equation_str = " + ".join(equation_parts)
full_symbolic_str = " + ".join(symbolic_parts)

# Print the final result, showing the formula and the calculated numbers.
print(f"The dimension of H_b^{n}(T x T; R) is calculated using the Künneth formula:")
print(f"dim H_b^{n}(T x T; R) = {full_symbolic_str}")
print(f"                   = {full_equation_str}")
print(f"                   = {total_dim}")
