import numpy as np

def calculate_relative_triplets(p):
    """
    Calculates a value proportional to the number of triplets required for a given alignment p.
    This is based on the inverse of the mutual information between the teacher and student.
    Note: For p=0 or p=1, the limit of p*log2(p) is 0. For p=0.5, the information is 0.
    """
    if p == 0.5:
        return float('inf')
    if p == 0 or p == 1:
        return 1.0
    
    # Mutual Information = 1 - BinaryEntropy(p)
    # BinaryEntropy(p) = -p*log2(p) - (1-p)*log2(1-p)
    information_gain = 1 - (-(p * np.log2(p) + (1 - p) * np.log2(1 - p)))
    
    return 1 / information_gain

def main():
    """
    Demonstrates the relationship between alignment (p) and the number of required triplets.
    """
    print("This script calculates the relative number of triplets required for learning, based on the probabilistic representational alignment 'p'.")
    print("The number of triplets is inversely proportional to the information gained from each triplet.")
    print("Number of Triplets ∝ 1 / (1 + p*log2(p) + (1-p)*log2(1-p))\n")

    # Case 1: High Anti-Alignment (p is close to 0)
    p_low = 0.01
    triplets_low = calculate_relative_triplets(p_low)
    print(f"Case 1: High Anti-Alignment (p = {p_low})")
    print("The student can reliably invert the teacher's statements, so information gain is high and few triplets are needed.")
    print(f"Num_Triplets ∝ 1 / (1 + {p_low}*log2({p_low}) + (1-{p_low})*log2(1-{p_low})) = {triplets_low:.4f}\n")

    # Case 2: Random Alignment (p = 0.5)
    p_mid = 0.5
    triplets_mid_str = "Infinity"
    print(f"Case 2: Random Alignment (p = {p_mid})")
    print("The teacher's statements are random noise to the student, providing no information. The number of triplets required is theoretically infinite.")
    print(f"Num_Triplets ∝ 1 / (1 + {p_mid}*log2({p_mid}) + (1-{p_mid})*log2(1-{p_mid})) = 1 / 0 => {triplets_mid_str}\n")
    
    # Case 3: High Alignment (p is close to 1)
    p_high = 0.99
    triplets_high = calculate_relative_triplets(p_high)
    print(f"Case 3: High Alignment (p = {p_high})")
    print("The student can reliably trust the teacher's statements, so information gain is high and few triplets are needed.")
    print(f"Num_Triplets ∝ 1 / (1 + {p_high}*log2({p_high}) + (1-{p_high})*log2(1-{p_high})) = {triplets_high:.4f}\n")

    print("Conclusion: The number of triplets is low at the extremes (p=0 and p=1) and maximal in the middle (p=0.5). This describes a convex U-shaped relationship.")

if __name__ == "__main__":
    main()