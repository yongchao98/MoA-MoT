import numpy as np

def run_analysis():
    """
    This function analyzes and demonstrates the relationship between probabilistic
    representational alignment (p) and the number of triplets required for a
    student to learn the location of a new object.

    The core idea is that the number of triplets needed is inversely proportional
    to the information gain from each triplet.
    """

    print("This script analyzes the relationship between probabilistic alignment (p)")
    print("and the number of triplets needed for learning.")
    print("\nPlan:")
    print("1. The information gain per triplet is calculated using the binary entropy function: Gain = 1 - H(p).")
    print("2. The relative number of triplets is calculated as: 1 / Gain.")
    print("3. We compute this value for a range of 'p' from 0 to 1 to see the shape of the relationship.")
    print("-" * 70)
    print("The final equation for the relative number of triplets N(p) is:")
    print("N(p) = 1 / (1 - [-p*log2(p) - (1-p)*log2(1-p)])")
    print("-" * 70)
    print(f"{'Alignment (p)':<20} | {'Relative Number of Triplets':<30}")
    print("-" * 70)

    # Generate a range of p values from 0 to 1
    p_values = np.linspace(0, 1, 21)

    for p in p_values:
        # At p=0.5, entropy is 1, information gain is 0, so N is infinite.
        if np.isclose(p, 0.5):
            num_triplets_str = "inf (Maximum)"
            print(f"{p:<20.2f} | {num_triplets_str:<30}")
            continue

        # At p=0 or p=1, entropy is 0. Using the limit definition.
        if np.isclose(p, 0.0) or np.isclose(p, 1.0):
            binary_entropy = 0.0
        else:
            # H(p) = -p*log2(p) - (1-p)*log2(1-p)
            binary_entropy = -p * np.log2(p) - (1 - p) * np.log2(1 - p)

        information_gain = 1.0 - binary_entropy

        # Relative number of triplets is inversely proportional to information gain
        relative_triplets = 1.0 / information_gain

        # Print each number in the final equation for this value of p
        print(f"For p={p:.2f}, Relative Triplets = 1.0 / (1.0 - {binary_entropy:.4f}) = {relative_triplets:.4f}")


    print("-" * 70)
    print("As the table shows, the number of triplets is high for p near 0.5 and")
    print("low for p near 0 and 1. This forms a convex U-shaped curve.")

# Run the analysis
run_analysis()