import math

# The problem asks for the value of av_n^k(1324) where n=333 and k=3.
# This represents the number of permutations of length 333 with 3 inversions
# that avoid the pattern 1324.

# For a fixed number of inversions k, this value is known to stabilize for large n.
# Since n=333 is significantly larger than k=3, we can use this stable value.
# The stable value for k inversions is given by the (k+1)-th Catalan number, C_{k+1}.

# We need to calculate the Catalan number for k=3.
k = 3
catalan_index = k + 1

# The formula for the m-th Catalan number (starting with m=0) is:
# C_m = (1 / (m + 1)) * C(2m, m), where C(a, b) is "a choose b".
m = catalan_index

# We calculate C_4.
two_m = 2 * m
binomial_coefficient = math.comb(two_m, m)
catalan_number = binomial_coefficient // (m + 1)

print("The value of av_333^3(1324) is given by the (3+1)-th Catalan number, C_4.")
print(f"The formula for C_m is (1/(m+1)) * (2m choose m).")
print(f"For m = {m}, the formula is:")
print(f"C_{m} = (1 / ({m} + 1)) * ({two_m} choose {m})")
print(f"First, we calculate the binomial coefficient ({two_m} choose {m}):")
print(f"({two_m} choose {m}) = {binomial_coefficient}")
print(f"Then, we divide by (m + 1) = {m + 1}:")
print(f"C_{m} = {binomial_coefficient} / {m + 1}")
print(f"The final result is: {catalan_number}")
