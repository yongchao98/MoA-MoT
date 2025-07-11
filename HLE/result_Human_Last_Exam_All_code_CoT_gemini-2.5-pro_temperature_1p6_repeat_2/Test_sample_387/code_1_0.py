import numpy as np

# Let P be the point with log structure given by the monoid M = N^3.
# The corresponding log scheme is (X, M_X) where X = Spec(k) is a single point.
# The dimension of the underlying scheme X is 0.
dim_X = 0

# We consider P over a base S = Spec(k) with the trivial log structure.
# The dimension of the base S is 0.
dim_S = 0

# The log structure of P is given by the monoid N^3. The characteristic
# monoid M_X_bar is N^3.
# The rank of a free monoid N^k is k. So, the rank of N^3 is 3.
# We can represent the generators of N^3 as vectors in Z^3.
monoid_generators = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
# The rank is the dimension of the vector space spanned by the generators,
# which is the rank of the matrix of generators.
rank_M_X_bar = np.linalg.matrix_rank(monoid_generators)

# The base S has the trivial log structure, meaning its characteristic
# monoid is {1} (or {0} in additive notation), which has rank 0.
rank_M_S_bar = 0

# The log dimension of P is given by the formula:
# logdim(P) = dim(X) - dim(S) + rank(M_X_bar) - rank(M_S_bar)
log_dimension = dim_X - dim_S + rank_M_X_bar - rank_M_S_bar

# A log blowup is a birational morphism in log geometry that preserves
# the log dimension of the scheme.
# Therefore, the dimension of the log blowup of P is equal to the log dimension of P.

print("The dimension of the log blowup is the log dimension of the log point P.")
print("The formula for the log dimension of P is:")
print("dim(P_underlying) - dim(Base_underlying) + rank(P_monoid) - rank(Base_monoid)")
print(f"The equation with the given values is:")
print(f"{dim_X} - {dim_S} + {rank_M_X_bar} - {rank_M_S_bar} = {log_dimension}")
