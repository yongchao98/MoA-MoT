from math import comb

# Step 1: Define the number of number fields of degree n (from 1 to 5)
# unramified outside the prime 2. These numbers are based on known results
# from algebraic number theory databases (e.g., the LMFDB).
N1 = 1  # The field is Q
N2 = 3  # The fields are Q(sqrt(-1)), Q(sqrt(2)), Q(sqrt(-2))
N3 = 0  # There are no such cubic fields
N4 = 7  # There are 7 such quartic fields
N5 = 0  # There are no such quintic fields

print("This calculation determines the number of isomorphism classes of degree 5 del Pezzo surfaces over Q")
print("with good reduction everywhere except possibly at the prime 2.")
print("This is equivalent to counting 5-dimensional étale algebras over Q unramified outside {2}.\n")
print("The count is based on the number of number fields K/Q of degree n unramified outside {2}:")
print(f"Number of fields of degree 1 (N1): {N1}")
print(f"Number of fields of degree 2 (N2): {N2}")
print(f"Number of fields of degree 3 (N3): {N3}")
print(f"Number of fields of degree 4 (N4): {N4}")
print(f"Number of fields of degree 5 (N5): {N5}")
print("-" * 40)
print("We sum the number of possible algebras for each partition of 5:\n")

# Step 2: Calculate the number of algebras for each partition of 5.
# An algebra's structure corresponds to a partition of 5, e.g., a partition
# 4+1 corresponds to an algebra K_4 x K_1 where K_n is a field of degree n.

# Partition 5: A single quintic field. The number of such algebras is N5.
count_5 = N5
print(f"Partition 5: Corresponds to a quintic field. Count = N5 = {count_5}")

# Partition 4+1: A quartic field and a linear field. Count is N4 * N1.
count_4_1 = N4 * N1
print(f"Partition 4+1: Corresponds to K_4 x K_1. Count = N4 * N1 = {N4} * {N1} = {count_4_1}")

# Partition 3+2: A cubic field and a quadratic field. Count is N3 * N2.
count_3_2 = N3 * N2
print(f"Partition 3+2: Corresponds to K_3 x K_2. Count = N3 * N2 = {N3} * {N2} = {count_3_2}")

# Partition 3+1+1: A cubic field (the K_1 fields are just Q). Count is N3.
count_3_1_1 = N3
print(f"Partition 3+1+1: Corresponds to K_3 x Q x Q. Count = N3 = {count_3_1_1}")

# Partition 2+2+1: Two distinct quadratic fields. We choose 2 from N2 options.
# The algebra is K_2a x K_2b x Q. Count is C(N2, 2).
count_2_2_1 = comb(N2, 2)
print(f"Partition 2+2+1: Corresponds to K_2a x K_2b x Q. Count = C(N2, 2) = C({N2}, 2) = {count_2_2_1}")

# Partition 2+1+1+1: One quadratic field. The algebra is K_2 x Q x Q x Q.
# Count is N2.
count_2_1_1_1 = N2
print(f"Partition 2+1+1+1: Corresponds to K_2 x Q x Q x Q. Count = N2 = {count_2_1_1_1}")

# Partition 1+1+1+1+1: The trivial algebra Q x Q x Q x Q x Q. There is only one.
count_1_1_1_1_1 = 1
print(f"Partition 1+1+1+1+1: Corresponds to Q^5. Count = {count_1_1_1_1_1}")

print("-" * 40)

# Step 3: Sum the results to get the total count.
total_count = count_5 + count_4_1 + count_3_2 + count_3_1_1 + count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1

print("The total number of isomorphism classes is the sum of these counts.")
print(f"Final Equation: {count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1} = {total_count}")