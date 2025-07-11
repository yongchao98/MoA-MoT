import numpy as np

def calculate_information_cost():
    """
    Calculates and prints the relationship between probabilistic representational alignment (p)
    and the number of triplets needed for learning, based on information theory.
    """
    
    # The final equation relates the number of triplets (N) to the alignment probability (p).
    # N is inversely proportional to the channel capacity C(p).
    # N ∝ 1 / C(p)
    # The capacity C(p) is given by 1 minus the binary entropy H_b(p).
    # C(p) = 1 - H_b(p)
    # The binary entropy H_b(p) is calculated as:
    # H_b(p) = -p*log2(p) - (1-p)*log2(1-p)
    
    print("The relationship between the number of triplets needed (N) and alignment (p) is derived from information theory.")
    print("The key equation components are:")
    print("1. Binary Entropy: H_b(p) = -p * log2(p) - (1-p) * log2(1-p)")
    print("2. Channel Capacity: C(p) = 1 - H_b(p)")
    print("3. Number of Triplets: N(p) is proportional to 1 / C(p)\n")

    def binary_entropy(p):
        # Use an epsilon to avoid log(0) for p=0 or p=1
        epsilon = 1e-12
        p_safe = np.clip(p, epsilon, 1 - epsilon)
        return -p_safe * np.log2(p_safe) - (1 - p_safe) * np.log2(1 - p_safe)

    def channel_capacity(p):
        return 1 - binary_entropy(p)

    # p values to demonstrate the relationship
    p_values = [0.01, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.99]
    
    print("Below are the calculated values that show the shape of the relationship.")
    print("N(p) goes to infinity as p approaches 0.5.")
    print("-" * 55)
    print(f"{'p (Alignment)':<15} | {'H_b(p)':<10} | {'C(p) (Capacity)':<15} | {'N(p) (proportional)'}")
    print("-" * 55)
    
    for p in p_values:
        h_b = binary_entropy(p)
        # Add epsilon to capacity to avoid division by zero at p=0.5
        c_p = channel_capacity(p) + 1e-12 
        n_p = 1 / c_p
        
        # Outputting each number in the final equation calculation for this p
        # Final equation: N(p) = 1 / (1 - (-p*log2(p) - (1-p)*log2(1-p)))
        if p == 0.5:
            n_p_str = "Infinity"
        else:
            n_p_str = f"{n_p:.3f}"

        print(f"{p:<15.2f} | {h_b:<10.3f} | {c_p:<15.3f} | {n_p_str}")
    
    print("-" * 55)
    print("\nThe number of triplets is low at the extremes (p=0, p=1) and high in the middle (p=0.5).")
    print("This describes a convex U-shaped relationship.")


calculate_information_cost()