def count_divisors(n):
    """
    Calculates the number of positive divisors of a positive integer n.
    This is also known as the tau function, τ(n).
    """
    if n < 1:
        return 0
    count = 0
    for i in range(1, n + 1):
        if n % i == 0:
            count += 1
    return count

# Step 1: Identify the order of the Schur multiplier for the given simple group.
# For the simple group S = PSL(2, p) with prime p > 5, the Schur multiplier M(S)
# is the cyclic group of order 2.
order_of_schur_multiplier = 2

# Step 2: The total number of non-isomorphic covering groups of a simple group S
# is equal to the number of subgroups of its Schur multiplier M(S).

# Step 3: For a cyclic group of order n, the number of subgroups is equal to the
# number of divisors of n, τ(n).
# Since M(PSL(2,p)) is a cyclic group of order 2, we need to calculate τ(2).
num_coverings = count_divisors(order_of_schur_multiplier)

# Step 4: Print the reasoning and the final result.
print("The problem is to find the total number of smooth coverings of PSL(2, p) for p > 5.")
print("This corresponds to the number of non-isomorphic covering groups in group theory.")
print("This number is equal to the number of subgroups of the Schur multiplier of PSL(2, p).")
print(f"For p > 5, the Schur multiplier is a cyclic group of order {order_of_schur_multiplier}.")
print("The number of subgroups of a cyclic group of order n is τ(n), the number of divisors of n.")
print("\nFinal Calculation:")
print(f"τ({order_of_schur_multiplier}) = {num_coverings}")
print(f"\nThus, there are {num_coverings} such smooth coverings.")
