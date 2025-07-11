import math

def dim_O(k):
  """Calculates the dimension of the orthogonal group O(k)."""
  if k < 0:
    return 0
  return k * (k - 1) // 2

# Number of replicas
n = 2

# Calculate the dimensions of the component groups O(2n) and O(n)
dim_O_2n = dim_O(2 * n)
dim_O_n = dim_O(n)

# Calculate the dimension of the numerator and denominator bosonic subgroups
dim_numerator = 2 * dim_O_2n
dim_denominator = 4 * dim_O_n

# Calculate the final result
num_variables = dim_numerator - dim_denominator

# Print the calculation steps as an equation
print("The number of non-Grassman variables is the dimension of the bosonic manifold M_B.")
print("dim(M_B) = dim(O(2n) x O(2n)) - dim((O(n) x O(n)) x (O(n) x O(n)))")
print(f"For n = {n}:")
print(f"dim(M_B) = (2 * dim(O({2*n}))) - (4 * dim(O({n})))")
print(f"           = (2 * {dim_O_2n}) - (4 * {dim_O_n})")
print(f"           = {dim_numerator} - {dim_denominator}")
print(f"           = {num_variables}")
