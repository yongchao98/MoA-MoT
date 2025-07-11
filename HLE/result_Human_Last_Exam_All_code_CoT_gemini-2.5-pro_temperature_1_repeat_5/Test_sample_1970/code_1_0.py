import math

def check_finite_analogue():
    """
    Analyzes a finite analogue of the set-theoretic problem using the
    probabilistic method.
    """
    print("This script explores a finite analogue of the problem.")
    print("We want to see if a function f: {0..m-1}^2 -> {0..k-1} exists such that")
    print("for every subset X of size n, the image f''([X]^2) has size k.\n")

    try:
        m_str = input("Enter a finite value for m (analogous to kappa^{++}): ")
        m = int(m_str)
        n_str = input("Enter a finite value for n (size of subset x, analogous to kappa+): ")
        n = int(n_str)
        k_str = input("Enter a finite value for k (number of colors, analogous to kappa): ")
        k = int(k_str)

        if not (m > n > k > 1):
            print("Please ensure m > n > k > 1 for a non-trivial problem.")
            return

    except ValueError:
        print("Invalid input. Please enter integers.")
        return

    # Calculate components of the probability bound inequality:
    # C(m, n) * k * (1 - 1/k)**C(n, 2) < 1

    try:
        # Number of subsets of size n
        num_subsets = math.comb(m, n)

        # Number of pairs in a subset of size n
        num_pairs_in_subset = math.comb(n, 2)
    except ValueError:
        print(f"Cannot compute combinations for m={m}, n={n}. Values might be too large or invalid.")
        return

    # Probability that one specific color is missing from the image
    # on a specific subset
    prob_one_color_missing = (1 - 1/k) ** num_pairs_in_subset

    # Upper bound on the probability that the image is not full size (union bound)
    prob_fail_one_subset = k * prob_one_color_missing

    # Upper bound on the probability that there exists any failing subset
    total_fail_prob_bound = num_subsets * prob_fail_one_subset

    # --- Output the results and the final equation ---
    print("\n--- Analysis using the Probabilistic Method ---")
    print("We check if the following inequality holds:")
    print("P(failure) <= C(m, n) * k * (1 - 1/k)^C(n, 2) < 1\n")
    
    print("Substituting the given numbers:")
    print(f"m = {m}")
    print(f"n = {n}")
    print(f"k = {k}\n")
    
    print("The numbers in the final equation are:")
    # Using scientific notation for potentially large or small numbers
    print(f"Number of subsets C({m}, {n}) = {num_subsets:e}")
    print(f"Number of colors k = {k}")
    print(f"Base of exponent (1 - 1/k) = {1 - 1/k}")
    print(f"Exponent C({n}, {2}) = {num_pairs_in_subset}")
    
    print("\nFinal Probability Bound Calculation:")
    print(f"P(failure) <= {num_subsets:e} * {k} * ({1-1/k})^{num_pairs_in_subset}")
    print(f"P(failure) <= {total_fail_prob_bound:e}\n")

    if total_fail_prob_bound < 1:
        print("Result: The probability bound is less than 1.")
        print("This proves that for these finite parameters, there MUST exist at least one")
        print("function f with the desired property.")
    else:
        print("Result: The probability bound is not less than 1.")
        print("The probabilistic method is inconclusive for these parameters. A function")
        print("might still exist, but this method cannot prove it.")

if __name__ == '__main__':
    check_finite_analogue()