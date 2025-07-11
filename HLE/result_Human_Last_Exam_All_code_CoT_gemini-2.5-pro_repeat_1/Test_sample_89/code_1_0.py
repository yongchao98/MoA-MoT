import numpy as np
import matplotlib.pyplot as plt

def number_of_triplets_vs_alignment(p):
    """
    Calculates a quantity proportional to the number of triplets needed
    for a given alignment probability p.

    Args:
      p: A numpy array of probability values between 0 and 1.

    Returns:
      A numpy array representing the relative number of triplets.
    """
    # To avoid log(0) and division by zero, we clip p to be slightly away from 0, 1, and 0.5.
    p_safe = np.clip(p, 1e-9, 1 - 1e-9)

    # H(p) is the binary entropy function.
    h_p = -p_safe * np.log2(p_safe) - (1 - p_safe) * np.log2(1 - p_safe)
    
    # I(p) is the mutual information.
    i_p = 1 - h_p
    
    # The number of triplets N(p) is inversely proportional to the information I(p).
    # We add a small epsilon to I(p) to prevent division by zero at p=0.5.
    n_p = 1 / (i_p + 1e-9)
    
    return n_p

def main():
    """
    Main function to calculate, plot, and explain the relationship.
    """
    # Create a range of p values from 0 to 1.
    p_values = np.linspace(0, 1, 500)
    
    # Calculate the corresponding number of triplets.
    n_values = number_of_triplets_vs_alignment(p_values)

    print("Analysis of the relationship:")
    print("Let N(p) be the number of triplets needed for an alignment p.")
    print("-" * 30)
    print("p=0.0 (Perfect anti-alignment): N(p) is minimal. Information per triplet is maximal.")
    print("p=0.5 (Random alignment): N(p) is maximal (infinite). Information per triplet is zero.")
    print("p=1.0 (Perfect alignment): N(p) is minimal. Information per triplet is maximal.")
    print("-" * 30)
    print("The overall shape is a 'U', with minima at p=0 and p=1, and a vertical asymptote at p=0.5.")
    print("The choice between 'Convex U-shaped' and 'Concave U-shaped' depends on the curvature of the arms of the U.")
    print("A detailed mathematical analysis shows the function has vertical tangents at its minima (p=0 and p=1), a characteristic of a shape formed by concave curves.")

    # Plotting the relationship
    plt.figure(figsize=(10, 6))
    plt.plot(p_values, n_values, lw=2)
    
    # Set a practical y-limit to see the shape near the minima
    plt.ylim(0, 40)
    
    plt.xlabel("p (Probabilistic Representational Alignment)", fontsize=12)
    plt.ylabel("N(p) (Number of Triplets Needed)", fontsize=12)
    plt.title("Relationship between Alignment and Number of Triplets", fontsize=14)
    plt.grid(True)
    
    # Annotate the key features
    plt.axvline(x=0.5, color='r', linestyle='--', label='p=0.5 (Infinite Triplets)')
    plt.text(0.01, 15, 'Minimal Triplets\n(Perfect Anti-alignment)', ha='left', color='blue')
    plt.text(0.99, 15, 'Minimal Triplets\n(Perfect Alignment)', ha='right', color='blue')
    plt.text(0.5, 30, 'Zero Information', ha='center', color='red', rotation=90)
    
    # Annotate the shape of the arms
    plt.text(0.2, 5, 'Concave Arm', ha='center', fontsize=11, style='italic')
    plt.text(0.8, 5, 'Concave Arm', ha='center', fontsize=11, style='italic')
    
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
