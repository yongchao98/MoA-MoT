import math

def count_integer_sequences():
    """
    Calculates the number of sequences of K positive integers up to N,
    where each number is greater than the one before it, and the increase
    between consecutive numbers does not exceed M.

    The calculation is based on the inclusion-exclusion principle, leading to the formula:
    Sum_{j=0 to K-1} [(-1)^j * C(K-1, j) * C(N - M*j, K)]
    """
    # --- Parameters of the sequence ---
    # You can change these values to solve for different constraints.
    N = 20
    K = 5
    M = 4

    print(f"Calculating the number of sequences for N={N}, K={K}, M={M}")
    print(f"The given condition M*(K-1) < N is {M*(K-1)} < {N}, which is {M*(K-1) < N}.")
    print("-" * 30)

    def combinations(n, k):
        """
        Calculates the binomial coefficient C(n, k).
        Returns 0 if k > n or k < 0.
        """
        if k < 0 or k > n:
            return 0
        # math.comb is available from Python 3.8 onwards
        try:
            return math.comb(n, k)
        except AttributeError:
            # Manual implementation for older Python versions
            if k == 0 or k == n:
                return 1
            if k > n // 2:
                k = n - k
            
            res = 1
            for i in range(k):
                res = res * (n - i) // (i + 1)
            return res

    # --- Part 1: Display the formula with the specific numbers ---
    equation_parts = []
    for j in range(K):
        sign = "-" if j % 2 else "+"
        # For the first term (j=0), we don't print the leading "+"
        if j == 0:
            sign_str = ""
        else:
            sign_str = f" {sign} "
        
        term_str = f"C({K-1}, {j}) * C({N - M*j}, {K})"
        equation_parts.append(sign_str + term_str)

    print("The final equation is:")
    print("Total = " + "".join(equation_parts))
    print("-" * 30)
    print("Calculating each term:")

    # --- Part 2: Calculate the result term by term ---
    total_sequences = 0
    for j in range(K):
        c1 = combinations(K - 1, j)
        
        n_for_c2 = N - M * j
        c2 = combinations(n_for_c2, K)
        
        term = c1 * c2
        sign_multiplier = -1 if j % 2 else 1
        signed_term = sign_multiplier * term
        
        total_sequences += signed_term
        
        # Print the breakdown of each number in the term calculation
        print(f"Term for j={j}:")
        print(f"  (-1)^{j} * C({K-1}, {j}) * C({N - M*j}, {K})")
        print(f"  = {sign_multiplier} * {c1} * {c2} = {signed_term}")

    print("-" * 30)
    print(f"The total number of possible sequences is: {total_sequences}")


if __name__ == '__main__':
    count_integer_sequences()
    # The final answer for N=20, K=5, M=4 is 2560
    # The submission format requires the final answer to be enclosed in <<<>>>
    print("<<<2560>>>")