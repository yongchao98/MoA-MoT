# This script computes the dimension of the degree 4 bounded cohomology group
# of the direct product of Thompson's group T with itself, with trivial real coefficients.

# The computation is based on the Künneth formula for bounded cohomology and
# known results about the bounded cohomology of Thompson's group T.

# Let d_k be the dimension of the k-th bounded cohomology group of T with real coefficients.
# The known dimensions are:
# d_0 = 1
# d_1 = 0
# d_2 = 1
# d_k = 0 for k >= 3
dims = {
    0: 1,
    1: 0,
    2: 1,
    3: 0,
    4: 0,
}

# We want to compute the dimension for degree n=4.
n = 4

# According to the Künneth formula, the dimension is the sum of products:
# dim H_b^4(T x T) = d_0*d_4 + d_1*d_3 + d_2*d_2 + d_3*d_1 + d_4*d_0

terms = []
total_dimension = 0

for p in range(n + 1):
    q = n - p
    dim_p = dims.get(p, 0)
    dim_q = dims.get(q, 0)
    
    # Calculate the product for the current term in the sum
    term_value = dim_p * dim_q
    
    # Add the value to the total dimension
    total_dimension += term_value
    
    # Store the string representation of the term for the final equation
    terms.append(f"{dim_p} * {dim_q}")

# Construct the full equation string
equation_string = " + ".join(terms)

# Print the final equation with each number and the result
print(f"The dimension is computed by the following sum:")
print(f"{equation_string} = {total_dimension}")