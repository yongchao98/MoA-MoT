def count_divisors(n):
    """
    Calculates the number of positive divisors for a given integer n.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    
    count = 0
    divisors_list = []
    # A number's divisors can be found by checking all integers from 1 to n.
    for i in range(1, n + 1):
        if n % i == 0:
            count += 1
            divisors_list.append(i)
    return count, divisors_list

def solve_coverings_problem():
    """
    Solves the problem by calculating the number of subgroups of the Schur multiplier.
    """
    # For the simple group S = PSL(2, p) with p > 5 prime, the Schur multiplier
    # M(S) is the cyclic group of order 2, denoted Z_2.
    n = 2

    # The number of non-isomorphic covering groups of S is the number of subgroups of M(S).
    # For a cyclic group Z_n, this is the number of divisors of n, tau(n).
    num_coverings, divisors = count_divisors(n)
    
    print(f"The number of smooth coverings corresponds to the number of subgroups of the Schur multiplier M(PSL(2, p)).")
    print(f"For p > 5, M(PSL(2, p)) is the cyclic group Z_{n}.")
    print(f"The number of subgroups of Z_{n} equals the number of divisors of {n}, which are: {', '.join(map(str, divisors))}.")
    print(f"The final calculation is:")
    print(f"Total number of coverings = tau({n}) = {num_coverings}")

if __name__ == "__main__":
    solve_coverings_problem()