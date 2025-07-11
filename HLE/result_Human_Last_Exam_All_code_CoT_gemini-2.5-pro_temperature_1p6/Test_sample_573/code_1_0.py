import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the given problem.

    The problem asks for the number of n-simplices in the over-category
    N_.(Z_N)_{k/}, which corresponds to counting non-decreasing sequences
    of n+1 integers in the range [-N, k].

    The formula for this is C(k + N + n + 1, n + 1), where C is the
    binomial coefficient ("choose" function).
    """
    N = 200
    k = 13

    print(f"Calculating the number of n-simplices for N={N}, k={k}, and n <= 5.")
    print("-" * 20)

    results = []
    for n in range(6):
        # The number of available integers is k - (-N) + 1 = k + N + 1
        # We are choosing n+1 integers with replacement.
        # The formula is Binomial( (k+N+1) + (n+1) - 1, n+1)
        top = k + N + n + 1
        bottom = n + 1
        
        # Using math.comb for combinations
        num_simplices = math.comb(top, bottom)
        results.append(num_simplices)

        print(f"For n = {n}, the number of simplices is given by the equation:")
        # Show the formula with the numbers plugged in
        equation_str = f"  C({k} + {N} + {n} + 1, {n} + 1)"
        # Show the evaluated numbers in the binomial coefficient
        evaluated_str = f"C({top}, {bottom})"
        # Show the final result
        result_str = f"{num_simplices}"
        
        print(f"{equation_str} = {evaluated_str} = {result_str}")

if __name__ == '__main__':
    solve_simplices_count()