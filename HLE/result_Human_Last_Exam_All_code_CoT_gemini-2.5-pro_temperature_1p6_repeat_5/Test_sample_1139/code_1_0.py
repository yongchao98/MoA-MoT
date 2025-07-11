# This script calculates the number of non-Grassman variables for a specific supersymmetric sigma-model.

# The problem specifies n=2 replicas for symmetry class D.
n = 2

# The bosonic manifold for this model is O(2n)/U(n).
# So we need to find the dimension of O(p)/U(q) where p=2n and q=n.
p = 2 * n
q = n

# The dimension of the orthogonal group O(p) is p*(p-1)/2.
dim_O_p = p * (p - 1) // 2

# The dimension of the unitary group U(q) is q^2.
dim_U_q = q**2

# The dimension of the manifold is the difference between the dimensions of the two groups.
result = dim_O_p - dim_U_q

# We print the final equation showing how the result is derived.
print(f"The number of variables is given by the dimension of the manifold O({p})/U({q}).")
print(f"This is calculated as dim(O({p})) - dim(U({q})).")
print(f"The final calculation is: {dim_O_p} - {dim_U_q} = {result}")
