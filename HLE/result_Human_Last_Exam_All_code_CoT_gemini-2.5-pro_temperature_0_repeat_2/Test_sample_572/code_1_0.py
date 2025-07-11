import math

def H2(delta):
    """
    Computes the binary entropy function H_2(delta).
    This is used in information theory and combinatorics to bound the size of sets.
    H_2(delta) = -delta * log2(delta) - (1 - delta) * log2(1 - delta)
    """
    # Handle edge cases where the formula is undefined but the limit is 0.
    if delta == 0 or delta == 1:
        return 0
    if delta < 0 or delta > 1:
        raise ValueError("Delta must be a probability, i.e., between 0 and 1.")
    
    return -delta * math.log2(delta) - (1 - delta) * math.log2(1 - delta)

def calculate_max_rank_factor(delta):
    """
    Calculates the constant factor 'c' for the rank r = c*N.
    This is derived from the probabilistic proof that random matrices are rigid.
    A matrix is likely non-rigid if the number of low-rank matrices plus their
    "neighborhoods" covers the entire space of matrices. This gives the inequality:
    H_2(delta) < (1 - c)^2, where r = c*N.
    """
    # The problem assumes delta is a small constant.
    # If delta >= 0.5, the concept of rigidity is less meaningful as one can
    # change half the entries to make the matrix zero (rank 0).
    if delta <= 0 or delta >= 0.5:
        return 0
    
    h = H2(delta)
    # From H_2(delta) < (1-c)^2, we get sqrt(H_2(delta)) < 1-c, so c < 1 - sqrt(H_2(delta)).
    c = 1 - math.sqrt(h)
    return c

def main():
    """
    Main function to explain and demonstrate the calculation.
    """
    print("Step 1: Understanding the Algorithm's Power")
    print("An FNP algorithm with an NP oracle can solve any search problem in NP.")
    print("The problem of finding a rigid matrix is a search problem for a co-NP set, as verifying rigidity is in co-NP.")
    print("Assuming a powerful enough oracle (e.g., if NP=co-NP, or a slightly stronger oracle), we can construct such a matrix if one is known to exist.")
    print("\nStep 2: Finding the Existence Limit")
    print("The question thus becomes: what is the largest 'r' for which a (delta, r)-rigid matrix is known to exist?")
    print("The best existence bounds come from the probabilistic method.")
    
    # The problem states delta is a small constant. Let's use a sample value.
    delta = 0.01
    
    print(f"\nStep 3: Calculating the Bound for delta = {delta}")
    c = calculate_max_rank_factor(delta)
    h_val = H2(delta)
    
    print("The bound is derived from the inequality involving the binary entropy function H_2(delta).")
    print(f"For delta = {delta}, the binary entropy H_2({delta}) is approximately {h_val:.4f}.")
    
    print("\nThe largest rank 'r' is of the form r = c * N, where 'c' must satisfy c < 1 - sqrt(H_2(delta)).")
    print("The final equation for the maximum factor 'c' is:")
    print(f"c < 1 - sqrt({h_val:.4f})")
    print(f"c < {1 - math.sqrt(h_val):.4f}")
    
    print("\nConclusion:")
    print("The largest rank 'r' for which such a matrix can be constructed is a linear function of N.")
    print(f"For delta = {delta}, the rank can be up to r ≈ {c:.4f} * N.")

if __name__ == "__main__":
    main()