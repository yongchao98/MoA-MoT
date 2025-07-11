import math

def count_divisors(n):
    """
    Calculates the number of divisors of a positive integer n.
    This is also known as the tau function, τ(n).
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("Input must be a positive integer.")
    
    count = 0
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            # If divisors are different, add two
            if n // i == i:
                count += 1
            else:
                count += 2
    return count

# For the simple group S = PSL(2, p) with p > 5 prime,
# the Schur multiplier M(S) is the cyclic group Z_2.
# The order of the Schur multiplier is 2.
order_of_schur_multiplier = 2

# The total number of non-isomorphic quasi-simple covering groups of S
# is equal to the number of subgroups of its Schur multiplier M(S).
# For a cyclic group Z_n, the number of subgroups is τ(n),
# the number of divisors of n.
total_coverings = count_divisors(order_of_schur_multiplier)

# Print the reasoning and the final equation.
# The equation is: Number of Coverings = τ(|M(S)|) = τ(2)
print("The number of quasi-simple covering groups of PSL(2, p) is the number of subgroups of its Schur Multiplier.")
print("For p > 5 prime, the Schur Multiplier of PSL(2, p) is Z_2, which has order 2.")
print(f"The number of subgroups of Z_2 is τ({order_of_schur_multiplier}), which is the number of divisors of {order_of_schur_multiplier}.")
print(f"Final Equation: Total Coverings = τ({order_of_schur_multiplier}) = {total_coverings}")
