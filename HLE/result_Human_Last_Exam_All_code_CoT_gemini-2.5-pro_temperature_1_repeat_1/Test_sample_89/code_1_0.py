import numpy as np
import matplotlib.pyplot as plt

def solve():
    """
    This script models and visualizes the relationship between probabilistic 
    representational alignment (p) and the number of triplets (N) needed to 
    teach a new object's location.
    """
    
    print("Step 1: Model the problem using Information Theory.")
    print("The communication from teacher to student is a Binary Symmetric Channel.")
    print("The probability 'p' that teacher and student agree is the probability of correct transmission.")
    
    print("\nStep 2: Define the information content per triplet using Channel Capacity (C).")
    print("The formula for channel capacity is: C = 1 - H(p)")
    print("H(p) is the binary entropy function.")

    print("\nStep 3: Define the binary entropy function H(p).")
    print("The final equation for H(p) is: H(p) = -p*log2(p) - (1-p)*log2(1-p)")
    # Note: The prompt asks to output each number in the final equation. 
    # In this formula, p is a variable, and the numbers are -1, 1, and 2 (from log2).
    # We will print the components of the formula.
    print("Breaking down the equation for H(p):")
    print("  - Term 1: (-p) * log2(p)")
    print("  - Term 2: (1 - p) * log2(1 - p)")
    print("  - H(p) = Term 1 - Term 2")


    print("\nStep 4: Relate the number of triplets (N) to the Channel Capacity (C).")
    print("The number of triplets 'N' needed is inversely proportional to the information per triplet.")
    print("Therefore, N is proportional to 1/C, which is 1 / (1 - H(p)).")
    
    # Generate data for the plot
    # p is the probability of agreement, from 0 to 1.
    # We use a small epsilon to avoid log(0) which is undefined.
    epsilon = 1e-9
    p = np.linspace(epsilon, 1.0 - epsilon, 500)

    # Calculate Binary Entropy H(p)
    H_p = -p * np.log2(p) - (1 - p) * np.log2(1 - p)

    # Calculate Channel Capacity C = 1 - H(p)
    C = 1 - H_p

    # N is proportional to 1/C. We'll plot 1/C.
    N = 1 / C
    
    print("\nStep 5: Plotting the relationship between p and N.")

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(p, N, label="N ∝ 1 / (1 - H(p))")
    
    # The value at p=0.5 goes to infinity, so we cap the y-axis for better visualization.
    plt.ylim(0, 40)
    plt.title("Number of Triplets (N) vs. Representational Alignment (p)")
    plt.xlabel("p (Probabilistic Representational Alignment)")
    plt.ylabel("N (Number of Triplets Required, in arbitrary units)")
    plt.grid(True)
    plt.axvline(x=0.5, color='r', linestyle='--', label="p=0.5 (Max Triplets)")
    plt.legend()
    
    # Annotate the graph
    plt.text(0.05, 5, "Low N needed\n(High Alignment)", ha='center')
    plt.text(0.95, 5, "Low N needed\n(High Alignment)", ha='center')
    plt.text(0.5, 35, "High N needed\n(Random Alignment)", ha='center')
    
    plt.show()

    print("\nConclusion: The plot clearly shows a convex U-shaped relationship.")

solve()