def count_divisors(n):
    """
    Calculates the number of positive divisors of an integer n.
    Returns the count and the list of divisors.
    """
    if n <= 0:
        return 0, []
    
    divisor_count = 0
    divisors_list = []
    for i in range(1, n + 1):
        if n % i == 0:
            divisor_count += 1
            divisors_list.append(i)
            
    return divisor_count, divisors_list

def solve_covering_groups_problem():
    """
    Solves the problem by finding the number of covering groups of PSL(2, p).
    """
    # The problem is to find the total number of non-isomorphic quasi-simple groups G
    # such that G/Z(G) is isomorphic to the simple group S = PSL(2, p) for a prime p > 5.
    # This number is equal to the number of subgroups of the Schur multiplier of S, M(S).
    
    # For S = PSL(2, p) with p > 5, the Schur multiplier M(S) is a cyclic group of order 2.
    order_of_schur_multiplier = 2
    
    # The number of subgroups of a cyclic group of order n is the number of divisors of n.
    # We need to calculate the number of divisors of 2.
    num_subgroups, divisors = count_divisors(order_of_schur_multiplier)
    
    print("Step 1: The total number of smooth coverings corresponds to the number of non-isomorphic covering groups of PSL(2, p).")
    print("Step 2: This number is equal to the number of subgroups of the Schur multiplier, M(PSL(2, p)).")
    print(f"Step 3: For a prime p > 5, the order of the Schur multiplier is {order_of_schur_multiplier}.")
    print(f"Step 4: The number of subgroups is the number of divisors of {order_of_schur_multiplier}.")
    
    # As requested, showing the numbers in the final equation/calculation.
    # The equation is: Number of coverings = Number of divisors of 2.
    print(f"The divisors of {order_of_schur_multiplier} are {divisors}.")
    print(f"The number of divisors of {order_of_schur_multiplier} is {num_subgroups}.")
    
    print("\nConclusion:")
    print(f"The total number of such smooth coverings is {num_subgroups}.")

solve_covering_groups_problem()