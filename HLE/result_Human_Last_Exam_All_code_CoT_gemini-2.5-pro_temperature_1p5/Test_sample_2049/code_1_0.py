import math

# Step 1: Explain the problem in number theoretic terms.
# The problem is equivalent to counting the number of degree-5 étale algebras over Q
# which are unramified outside the set of primes S = {2, infinity}.
# This can be calculated by summing the counts for each partition of 5.

# Step 2: Use known data for the number of relevant number fields.
# Let a_n be the number of number fields of degree n unramified outside {2, infinity}.
# These values are based on the Jones-Roberts database.
a1 = 1  # The rational field Q
a2 = 3  # The quadratic fields Q(sqrt(-1)), Q(sqrt(2)), Q(sqrt(-2))
a3 = 2  # Two non-isomorphic cubic fields
a4 = 9  # Nine non-isomorphic quartic fields
a5 = 10 # Ten non-isomorphic quintic fields

# Step 3: Calculate the number of algebras for each partition of 5.

# Partition 5: Corresponds to a single quintic field.
# Number of choices = a5
n_5 = a5

# Partition 4+1: Corresponds to a product of a quartic field and Q.
# Number of choices = a4 * a1
n_41 = a4 * a1

# Partition 3+2: Corresponds to a product of a cubic and a quadratic field.
# Number of choices = a3 * a2
n_32 = a3 * a2

# Partition 3+1+1: Corresponds to a product of a cubic field and Q x Q.
# Number of choices = a3
n_311 = a3

# Partition 2+2+1: Product of two quadratic fields and Q. The two fields are
# chosen from the set of a2=3 fields, with replacement.
# This is a combination with replacement problem: C(n+k-1, k) for n=a2, k=2.
# C(3+2-1, 2) = C(4, 2) = 6.
n_221 = math.comb(a2 + 2 - 1, 2)

# Partition 2+1+1+1: Product of one quadratic field and Q x Q x Q.
# Number of choices = a2
n_2111 = a2

# Partition 1+1+1+1+1: The trivial algebra Q x Q x Q x Q x Q.
# There is only one such algebra.
n_11111 = 1

# Step 4: Sum the counts for all partitions.
total = n_5 + n_41 + n_32 + n_311 + n_221 + n_2111 + n_11111

# Step 5: Print the final equation.
# The problem asks for the sum to be spelled out.
print(f"The total number is the sum of counts for each partition of 5:")
print(f"5: {n_5}")
print(f"4+1: {n_41}")
print(f"3+2: {n_32}")
print(f"3+1+1: {n_311}")
print(f"2+2+1: {n_221}")
print(f"2+1+1+1: {n_2111}")
print(f"1+1+1+1+1: {n_11111}")
print("Final Equation:")
print(f"{n_5} + {n_41} + {n_32} + {n_311} + {n_221} + {n_2111} + {n_11111} = {total}")
