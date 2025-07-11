import numpy as np
import warnings

def analyze_triplet_efficiency():
    """
    This script models the relationship between representational alignment (p)
    and the number of triplets required to teach a new concept.

    The analysis is based on information theory, where the number of triplets
    needed is inversely proportional to the information content of each triplet.
    """

    # We assume a fixed amount of total information is needed to learn the location.
    # The value itself doesn't affect the shape of the relationship, only the scale.
    TOTAL_INFO_NEEDED = 100.0  # bits

    def binary_entropy(p):
        """Calculates the binary entropy H(p) = -p*log2(p) - (1-p)*log2(1-p)."""
        # Handle edge cases p=0 and p=1 to avoid log(0) which is -inf.
        # The limit of p*log2(p) as p->0 is 0.
        if p == 0 or p == 1:
            return 0
        return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

    def number_of_triplets(p):
        """
        Calculates the number of triplets needed as a function of alignment p.
        This is total info needed / info per triplet.
        Info per triplet is 1 - H(p), the capacity of a binary symmetric channel.
        """
        # Information from one triplet (channel capacity)
        info_per_triplet = 1 - binary_entropy(p)

        # Avoid division by zero when p=0.5, where info_per_triplet is 0.
        if info_per_triplet == 0:
            return float('inf')
        
        return TOTAL_INFO_NEEDED / info_per_triplet

    print("This script calculates the theoretical number of triplets required to communicate")
    print(f"a new location, assuming {TOTAL_INFO_NEEDED} bits of information are needed.")
    print("The calculation is shown for various values of representational alignment 'p'.\n")

    print(f"{'p (Alignment)':<15} | {'Equation: N = I_total / (1 - H(p))':<40} | {'Number of Triplets (N)':<25}")
    print("-" * 85)

    p_values = [0.0, 0.01, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.99, 1.0]

    for p in p_values:
        # For p=0 or p=1, H(p)=0. N = 100 / (1-0) = 100
        # For p=0.5, H(p)=1. N = 100 / (1-1) = infinity
        required_triplets = number_of_triplets(p)
        
        info_gain = 1 - binary_entropy(p)
        equation_str = f"N = {TOTAL_INFO_NEEDED:.0f} / (1 - {binary_entropy(p):.3f}) = {TOTAL_INFO_NEEDED:.0f} / {info_gain:.3f}"

        print(f"{p:<15.2f} | {equation_str:<40} | {required_triplets:<25.2f}")

    print("\nObservation:")
    print("The number of triplets is minimal (100.00) at p=0.0 and p=1.0 (perfect alignment or anti-alignment).")
    print("The number increases as p moves towards 0.5.")
    print("At p=0.5 (random alignment), the number of triplets required is infinite.")
    print("This demonstrates a distinct U-shaped curve that is symmetric around p=0.5.")
    print("This relationship is a convex U-shape.")

if __name__ == '__main__':
    analyze_triplet_efficiency()
