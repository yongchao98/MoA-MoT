import math

def combinations_with_repetition(n, k):
    """
    Calculates the number of ways to choose k items from a set of n items with repetition.
    This is also known as "n choose k with repetition" or multiset coefficient.
    The formula is C(n + k - 1, k).
    """
    if k < 0 or n < 0:
        return 0
    if n == 0 and k > 0:
        return 0
    if k == 0:
        return 1
    return math.comb(n + k - 1, k)

# Step 1: Define the number of number fields unramified outside the prime 2 for each degree.
# These numbers are based on data from number field databases (e.g., LMFDB).
# n_d is the number of fields of degree d ramified only at 2.
# A key theorem states that such fields must have a degree that is a power of 2.
n_d = {
    1: 1,  # The field Q
    2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
    3: 0,  # Degree is not a power of 2
    4: 7,
    5: 0,  # Degree is not a power of 2
}

# Step 2: Calculate the number of algebras for each partition of 5.
# A quintic etale algebra is a product of number fields where the sum of degrees is 5.
# We sum the possibilities over all partitions of 5.

# Partition 5: A single quintic field.
# Number of choices = n_5
p5 = n_d[5]

# Partition 4+1: A quartic field and a linear field (Q).
# Number of choices = n_4 * n_1
p4_1 = n_d[4] * n_d[1]

# Partition 3+2: A cubic field and a quadratic field.
# Number of choices = n_3 * n_2
p3_2 = n_d[3] * n_d[2]

# Partition 3+1+1: A cubic field and two copies of Q.
# Corresponds to choosing 1 cubic field and a multiset of 2 linear fields.
p3_1_1 = n_d[3] * combinations_with_repetition(n_d[1], 2)

# Partition 2+2+1: Two quadratic fields and Q.
# Corresponds to choosing a multiset of 2 quadratic fields and 1 linear field.
p2_2_1 = combinations_with_repetition(n_d[2], 2) * n_d[1]

# Partition 2+1+1+1: One quadratic field and three copies of Q.
# Corresponds to choosing 1 quadratic field and a multiset of 3 linear fields.
p2_1_1_1 = n_d[2] * combinations_with_repetition(n_d[1], 3)

# Partition 1+1+1+1+1: Five copies of Q.
# Corresponds to choosing a multiset of 5 linear fields.
p1_1_1_1_1 = combinations_with_repetition(n_d[1], 5)

# Step 3: Sum the results and print the detailed calculation.
total = p5 + p4_1 + p3_2 + p3_1_1 + p2_2_1 + p2_1_1_1 + p1_1_1_1_1

print("The number of isomorphism classes is the sum of counts for each partition of 5:")
print(f"Partition 5: {p5}")
print(f"Partition 4+1: {p4_1}")
print(f"Partition 3+2: {p3_2}")
print(f"Partition 3+1+1: {p3_1_1}")
print(f"Partition 2+2+1: {p2_2_1}")
print(f"Partition 2+1+1+1: {p2_1_1_1}")
print(f"Partition 1+1+1+1+1: {p1_1_1_1_1}")
print(f"Total = {p5} + {p4_1} + {p3_2} + {p3_1_1} + {p2_2_1} + {p2_1_1_1} + {p1_1_1_1_1} = {total}")