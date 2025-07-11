import math

# We are looking for the number of internal adjunctions from [m] to [n].
# In this problem, m=23 and n=37.
m = 23
n = 37

# This problem is equivalent to counting the number of order-preserving maps
# R: [n] -> [m] such that R(n) = m.
# This in turn is equivalent to counting the number of non-decreasing sequences
# of length n (from R(0) to R(n-1)) with values in {0, ..., m}.

# The number of items to choose is k = n.
# The number of values to choose from is N = m + 1.
# By the stars and bars formula, the number of such sequences is C(k + N - 1, k).
# In our specific case, the condition R(n)=m reduces the problem to counting
# non-decreasing sequences of length n, from R(0) to R(n-1), with values in {0, ..., m}.
# Wait, the sequence is from R(0) to R(n-1), so its length is n.
# Let's re-state based on the derivation.
# We are counting maps from [n-1] to [m].
# k = n-1, l = m.
# Number of maps is C(l+k+1, k+1) = C(m + (n-1) + 1, (n-1)+1) = C(m+n, n).
# This gives C(23+37, 37) = C(60, 37)

# The number of items we are choosing is for R(0), ..., R(n-1), which is n values.
# The set of values they can take is {0, ..., m}, which has m+1 elements.
# k = n = 37
# N = m+1 = 24
# The number of combinations is C(k + N - 1, k) = C(37 + 24 - 1, 37) = C(60, 37).
# Note that C(60, 37) is the same as C(60, 60-37) = C(60, 23).

val_N = 60
val_k = 37

# Calculate the binomial coefficient C(60, 37)
result = math.comb(val_N, val_k)

# Print the final equation
print(f"The number of internal adjunctions is given by the binomial coefficient C({val_N}, {val_k}).")
print(f"C({val_N}, {val_k}) = {result}")
