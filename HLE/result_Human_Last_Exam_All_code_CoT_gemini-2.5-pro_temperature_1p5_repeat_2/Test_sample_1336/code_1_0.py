import math

def get_divisors(n):
    """Returns a sorted list of divisors of a positive integer n."""
    if n <= 0:
        return []
    divisors = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divisors.add(i)
            divisors.add(n // i)
    return sorted(list(divisors))

# Step 1: State the problem in terms of group theory.
# The number of "smooth coverings" corresponds to the number of quasi-simple covering groups of PSL(2,p).
# This number is equal to the number of subgroups of the Schur multiplier M(PSL(2,p)).

# Step 2: Identify the relevant Schur multiplier.
# For a prime p > 5, the Schur multiplier of PSL(2,p) is the cyclic group of order 2, C_2.
order_of_multiplier = 2

# Step 3: Relate the number of subgroups to the number of divisors.
# The number of subgroups of a cyclic group C_n is the number of divisors of n.
# So we need to find the number of divisors of 2.

# Step 4: Perform the calculation.
n = order_of_multiplier
divisors_list = get_divisors(n)
num_divisors = len(divisors_list)

# Step 5: Print the explanation and the result.
print(f"The total number of smooth coverings is equal to the number of subgroups of the Schur multiplier of PSL(2, p) for p > 5.")
print(f"The Schur multiplier for this group is the cyclic group C_n where n = {n}.")
print(f"The number of subgroups of C_{n} is equal to the number of divisors of n, which is denoted tau({n}).")
print("\nFinal Equation:")
print(f"The divisors of {n} are {divisors_list}.")
print(f"tau({n}) = Number of divisors = {num_divisors}")
print(f"\nTherefore, the total number of such smooth coverings is {num_divisors}.")
