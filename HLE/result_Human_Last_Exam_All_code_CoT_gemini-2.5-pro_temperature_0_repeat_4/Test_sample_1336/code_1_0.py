def count_divisors(n):
    """
    Calculates the number of positive divisors of a positive integer n.
    This is also the number of subgroups of the cyclic group Z_n.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("Input must be a positive integer.")
    
    count = 0
    divisors_list = []
    for i in range(1, n + 1):
        if n % i == 0:
            count += 1
            divisors_list.append(i)
    return count, divisors_list

def solve_covering_groups():
    """
    Solves the problem of finding the total number of smooth coverings for a given p.

    The logic is as follows:
    1. The total number of smooth coverings corresponds to the number of non-isomorphic
       covering groups of the simple group S = PSL(2, p).
    2. The number of covering groups of S is equal to the number of subgroups of its
       Schur Multiplier, M(S).
    3. For S = PSL(2, p) with p > 5, the Schur Multiplier is M(S) = Z_2 (the cyclic group of order 2).
    4. The number of subgroups of Z_n is equal to the number of divisors of n, tau(n).
    5. We need to calculate tau(2) to find the answer.
    """
    # The order of the Schur Multiplier M(PSL(2, p)) for p > 5 is 2.
    n = 2
    
    print("Step 1: The problem is to find the total number of smooth coverings.")
    print("This number is equivalent to the number of non-isomorphic covering groups of PSL(2, p).")
    print("-" * 20)
    
    print("Step 2: The number of covering groups of a simple group S is equal to the number of subgroups of its Schur multiplier, M(S).")
    print("-" * 20)

    print(f"Step 3: For S = PSL(2, p) with p > 5, the Schur multiplier is the cyclic group Z_{n}.")
    print("-" * 20)

    print(f"Step 4: The number of subgroups of Z_{n} is equal to the number of divisors of {n}, denoted as tau({n}).")
    print("-" * 20)
    
    print(f"Step 5: We calculate tau({n}).")
    num_subgroups, divisors = count_divisors(n)
    print(f"The divisors of {n} are: {divisors}.")
    
    # Output the final equation as requested
    print(f"The final equation is: tau({n}) = {num_subgroups}")
    print("-" * 20)
    
    print(f"Conclusion: The total number of such smooth coverings is {num_subgroups}.")

if __name__ == "__main__":
    solve_covering_groups()
