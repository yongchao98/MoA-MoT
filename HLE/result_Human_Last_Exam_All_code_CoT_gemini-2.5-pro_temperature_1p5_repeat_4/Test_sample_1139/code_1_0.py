# Number of replicas for the disordered system
n = 2

# The number of non-Grassman variables is the dimension of the bosonic manifold.
# For symmetry class D, this is the symmetric space G/H = O(2n)/(O(n)xO(n)).
# The dimension is calculated as dim(G) - dim(H).

# 1. Calculate the dimension of G = O(2n)
# The dimension of the orthogonal group O(k) is k*(k-1)/2.
k_G = 2 * n
dim_G = k_G * (k_G - 1) // 2

# 2. Calculate the dimension of H = O(n) x O(n)
# The dimension is dim(O(n)) + dim(O(n)).
k_H = n
dim_O_n = k_H * (k_H - 1) // 2
dim_H = 2 * dim_O_n

# 3. Calculate the final number of variables
num_variables = dim_G - dim_H

# Print the explanation and the final equation with all numbers
print(f"The number of bosonic variables for a class D system with n={n} replicas is:")
print(f"dim(O({k_G})) - dim(O({k_H}) x O({k_H}))\n")

print(f"First, calculate the dimension of O({k_G}):")
print(f"dim(O({k_G})) = {k_G} * ({k_G} - 1) / 2 = {dim_G}")

print(f"\nNext, calculate the dimension of O({k_H}) x O({k_H}):")
print(f"dim(O({k_H}) x O({k_H})) = 2 * (dim(O({k_H}))) = 2 * ({k_H} * ({k_H} - 1) / 2) = {dim_H}")

print("\nFinal Equation:")
# The final equation shows each number in the calculation
print(f"Number of variables = {dim_G} - {dim_H} = {num_variables}")