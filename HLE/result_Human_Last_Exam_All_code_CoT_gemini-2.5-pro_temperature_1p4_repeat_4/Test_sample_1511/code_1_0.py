import math

# Given parameters for the problem
m = 5
k = 2
t = 1 # for cross 1-intersecting families

# According to the Frankl-Tokushige theorem for cross-t-intersecting multisets,
# for m >= k+t, the maximal sum of sizes |F| + |G| is 2 * C(m+k-t-1, k-t).
# Here, we calculate this value.

# Calculate the arguments for the combinations function
numerator_arg = m + k - t - 1
denominator_arg = k - t

# Calculate the size of one of the maximal families
# This is C(5+2-1-1, 2-1) = C(5, 1)
size_of_one_family = math.comb(numerator_arg, denominator_arg)

# The maximal sum is twice this size
max_sum = 2 * size_of_one_family

print(f"The calculation for part (b) is based on the formula for the maximum sum of sizes of cross-intersecting families:")
print(f"Max Sum = 2 * C(m + k - t - 1, k - t)")
print(f"Substituting m={m}, k={k}, t={t}:")
print(f"Max Sum = 2 * C({m} + {k} - {t} - 1, {k} - {t}) = 2 * C({numerator_arg}, {denominator_arg})")
print(f"Max Sum = 2 * {size_of_one_family} = {max_sum}")