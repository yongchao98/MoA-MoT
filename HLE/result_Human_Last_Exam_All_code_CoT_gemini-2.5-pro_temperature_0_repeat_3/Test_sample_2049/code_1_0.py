import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k), also known as "n choose k".
    """
    if k < 0 or k > n:
        return 0
    # This is an efficient way to compute combinations without large intermediate numbers.
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

# Step 1: Explain the problem and method.
print("The number of isomorphism classes of del Pezzo surfaces of degree 5 over Q")
print("with good reduction everywhere except possibly at the prime 2 is equivalent to")
print("the number of isomorphism classes of etale algebras of degree 5 over Q")
print("that are unramified outside the set of primes S = {2, infinity}.\n")

# Step 2: Provide the necessary data.
# N_k is the number of number fields of degree k over Q unramified outside S.
# These values are obtained from number theory databases (e.g., LMFDB).
N = {
    1: 1,  # The field Q itself.
    2: 3,  # The fields Q(sqrt(-1)), Q(sqrt(2)), Q(sqrt(-2)).
    3: 0,  # There are no such cubic fields.
    4: 18, # From database look-up.
    5: 0   # There are no such quintic fields.
}

print("This number is calculated by summing contributions from each partition of 5.")
print("Let N_k be the number of number fields of degree k unramified outside {2, infinity}.")
print(f"The known values are: N_1={N[1]}, N_2={N[2]}, N_3={N[3]}, N_4={N[4]}, N_5={N[5]}.\n")

# Step 3: Calculate the number of etale algebras for each partition of 5.
# The number of etale algebras of degree n is given by summing over partitions of n.
# For a partition n = k1*m1 + k2*m2 + ..., the number is product of C(N_ki + mi - 1, mi).
total_count = 0
equation_parts = []

# Partition (5): Corresponds to field extensions of degree 5.
# Formula: N_5
count_5 = N[5]
total_count += count_5
equation_parts.append(str(count_5))
print(f"Contribution from partition (5): {count_5}")

# Partition (4, 1): Algebras of the form K_4 x K_1.
# Formula: N_4 * N_1
count_4_1 = N[4] * N[1]
total_count += count_4_1
equation_parts.append(str(count_4_1))
print(f"Contribution from partition (4,1): {count_4_1}")

# Partition (3, 2): Algebras of the form K_3 x K_2.
# Formula: N_3 * N_2
count_3_2 = N[3] * N[2]
total_count += count_3_2
equation_parts.append(str(count_3_2))
print(f"Contribution from partition (3,2): {count_3_2}")

# Partition (3, 1, 1): Algebras of the form K_3 x Q x Q.
# Formula: N_3 * C(N_1 + 2 - 1, 2)
count_3_1_1 = N[3] * combinations(N[1] + 1, 2)
total_count += count_3_1_1
equation_parts.append(str(count_3_1_1))
print(f"Contribution from partition (3,1,1): {count_3_1_1}")

# Partition (2, 2, 1): Algebras of the form K_2a x K_2b x Q.
# Formula: C(N_2 + 2 - 1, 2) * N_1
count_2_2_1 = combinations(N[2] + 1, 2) * N[1]
total_count += count_2_2_1
equation_parts.append(str(count_2_2_1))
print(f"Contribution from partition (2,2,1): {count_2_2_1}")

# Partition (2, 1, 1, 1): Algebras of the form K_2 x Q x Q x Q.
# Formula: N_2 * C(N_1 + 3 - 1, 3)
count_2_1_1_1 = N[2] * combinations(N[1] + 2, 3)
total_count += count_2_1_1_1
equation_parts.append(str(count_2_1_1_1))
print(f"Contribution from partition (2,1,1,1): {count_2_1_1_1}")

# Partition (1, 1, 1, 1, 1): The algebra Q x Q x Q x Q x Q.
# Formula: C(N_1 + 5 - 1, 5)
count_1_1_1_1_1 = combinations(N[1] + 4, 5)
total_count += count_1_1_1_1_1
equation_parts.append(str(count_1_1_1_1_1))
print(f"Contribution from partition (1,1,1,1,1): {count_1_1_1_1_1}")

# Step 4: Print the final sum and the result.
print("\nThe total number is the sum of these contributions.")
final_equation = " + ".join(equation_parts)
print(f"Total = {final_equation} = {total_count}")
