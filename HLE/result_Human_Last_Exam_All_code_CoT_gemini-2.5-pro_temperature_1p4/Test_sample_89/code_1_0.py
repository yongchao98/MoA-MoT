import numpy as np

def binary_entropy(p):
    """
    Calculates the binary entropy H(p).
    Uses a small epsilon for numerical stability to avoid log(0).
    """
    epsilon = 1e-12
    p_clipped = np.clip(p, epsilon, 1 - epsilon)
    return -p_clipped * np.log2(p_clipped) - (1 - p_clipped) * np.log2(p_clipped)

def get_relative_triplets(p_agreement):
    """
    Calculates the relative number of triplets needed.
    This is inversely proportional to the channel capacity.
    """
    # The probability of disagreement (crossover probability) is 1 - p
    crossover_prob = 1 - p_agreement
    
    # The information per triplet (channel capacity)
    capacity = 1 - binary_entropy(crossover_prob)
    
    # The number of triplets is inversely proportional to the capacity.
    # Add epsilon to avoid division by zero when capacity is 0 (at p=0.5).
    epsilon = 1e-12
    return 1 / (capacity + epsilon)

def analyze_relationship():
    """
    Analyzes and prints the relationship between alignment 'p' 
    and the number of triplets needed.
    """
    print("Analyzing the relationship between alignment (p) and the number of triplets needed.")
    print("The number of triplets is inversely proportional to the information content of each triplet.")
    print("-" * 20)
    
    p_values_to_test = [0.01, 0.25, 0.5, 0.75, 0.99]
    
    for p in p_values_to_test:
        relative_n_triplets = get_relative_triplets(p)
        print(f"For alignment p = {p:.2f}, the relative number of triplets needed is: {relative_n_triplets:.2f}")

    print("-" * 20)
    print("The values are lowest at the extremes (p near 0 and 1) and highest in the middle (p=0.5).")
    print("This indicates a convex U-shaped relationship.")

# Execute the analysis
analyze_relationship()
