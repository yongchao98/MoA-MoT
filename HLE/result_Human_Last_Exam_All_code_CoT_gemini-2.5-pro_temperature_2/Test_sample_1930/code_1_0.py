# The vector space of digitary functions can be decomposed into subspaces.

# Dimension of the subspace of constant functions, e.g., f(x) = c.
# This space is spanned by the function f(x) = 1.
dim_const = 1
print(f"Dimension of constant functions = {dim_const}")

# Dimension of the subspace of linear functions, f(x) = ax.
# This space is spanned by the function f(x) = x.
dim_linear = 1
print(f"Dimension of linear functions = {dim_linear}")

# Dimension of additional 'purely digitary' functions.
# These are related to functions h mapping digits D={0,...,9} to reals,
# with the constraint that h(0)=0 and h(9)=0.
# The dimension of the space of these h functions is the number of digits
# minus the two constraints.
# Digits are {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, so there are 10 digits.
num_digits = 10
num_constraints = 2
dim_h = num_digits - num_constraints
print(f"Dimension of h-based functions = {dim_h}")

# The total dimension is the sum of the dimensions of these subspaces.
total_dimension = dim_const + dim_linear + dim_h

print(f"\nThe total dimension is the sum of these parts:")
print(f"Total Dimension = {dim_const} + {dim_linear} + {dim_h} = {total_dimension}")
