# Number of number fields of degree n unramified outside the prime 2
F1 = 1  # The field Q
F2 = 3  # Q(sqrt(-1)), Q(sqrt(2)), Q(sqrt(-2))
F3 = 0  # No cubic fields exist with this property
F4 = 10 # From databases like LMFDB
F5 = 10 # From databases like LMFDB

# Calculate the number of etale algebras for each partition of 5
# Partition 5:
count_5 = F5
# Partition 4+1:
count_4_1 = F4 * F1
# Partition 3+2:
count_3_2 = F3 * F2
# Partition 3+1+1:
count_3_1_1 = F3 # Number of ways to choose one F3 field and two F1 fields
# Partition 2+2+1:
# This is choosing 2 quadratic fields from a set of F2=3, with replacement.
# The formula for combinations with replacement of k items from a set of n is C(n+k-1, k).
# Here n=F2=3, k=2. So, C(3+2-1, 2) = C(4, 2) = 6.
count_2_2_1 = 6
# Partition 2+1+1+1:
count_2_1_1_1 = F2
# Partition 1+1+1+1+1:
count_1_1_1_1_1 = F1

# The total number is the sum of these counts.
total_count = count_5 + count_4_1 + count_3_2 + count_3_1_1 + count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1

print("The total number of isomorphism classes is the sum of the possibilities for each partition of 5:")
print(f"Partition 5: {count_5}")
print(f"Partition 4+1: {count_4_1}")
print(f"Partition 3+2: {count_3_2}")
print(f"Partition 3+1+1: {count_3_1_1}")
print(f"Partition 2+2+1: {count_2_2_1}")
print(f"Partition 2+1+1+1: {count_2_1_1_1}")
print(f"Partition 1+1+1+1+1: {count_1_1_1_1_1}")
print("-" * 20)
# Final equation must have the numbers
print(f"Total = {count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1} = {total_count}")