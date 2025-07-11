import math

# This script calculates the number of nonzero terms in the asymptotic expansion
# by counting how many integers n from 1 to 100 satisfy the conditions for c_n != 0.

# Category 1: n = 2^p, where p is odd, and n <= 100.
p = 1
n = 2**p
cat1_nums = []
while n <= 100:
    cat1_nums.append(n)
    p += 2
    n = 2**p
count_cat1 = len(cat1_nums)

# Category 2: n = 2^p * m, where m is an odd integer > 1, p is even, and n <= 100.
# We sum the counts for p = 0, 2, 4, ...

# p = 0: n = m (odd integers from 3 to 99).
count_p0 = len(range(3, 100, 2))

# p = 2: n = 4 * m. We need 4m <= 100, so m <= 25.
# m must be an odd integer > 1.
# These are m = 3, 5, ..., 25.
count_p2 = len(range(3, 26, 2))

# p = 4: n = 16 * m. We need 16m <= 100, so m <= 6.25.
# m must be an odd integer > 1.
# These are m = 3, 5.
count_p4 = len(range(3, 6, 2))

# p = 6: n = 64 * m. We need 64m <= 100, so m <= 1.5625.
# There are no odd integers m > 1 that satisfy this. The count is 0.
# Higher even powers of p also result in a count of 0.

# The total number of non-zero terms is the sum of the counts from these categories.
total_nonzero_terms = count_cat1 + count_p0 + count_p2 + count_p4
counts = [count_cat1, count_p0, count_p2, count_p4]

# Print the final calculation as an equation, as requested.
equation_str = " + ".join(map(str, counts))
print(f"The number of non-zero terms is the sum of terms in each category:")
print(f"{equation_str} = {total_nonzero_terms}")
