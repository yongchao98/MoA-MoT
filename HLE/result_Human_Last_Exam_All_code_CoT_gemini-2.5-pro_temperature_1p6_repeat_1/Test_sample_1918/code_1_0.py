import math

def calculate_max_complementary_rank(d, n, r):
    """
    Calculates the maximal rank of the Choi matrix of the complementary channel.

    Args:
        d (int): Dimension of the input Hilbert space H1.
        n (int): Dimension of the output Hilbert space H2.
        r (int): Rank of the Choi matrix of the original channel Lambda.
    """
    # The existence of the channel Lambda imposes the constraint d <= n*r
    if d > n * r:
        print(f"Warning: No channel can exist for the given parameters as the constraint d <= n*r is not met.")
        print(f"d={d}, n={n}, r={r}, n*r={n*r}")
        print("Proceeding with calculation, but the premise is invalid.")
        
    print(f"Given parameters:")
    print(f"Dimension of input space (d): {d}")
    print(f"Dimension of output space (n): {n}")
    print(f"Rank of original Choi matrix (r): {r}")
    print("-" * 20)
    
    # The maximal rank of the complementary Choi matrix is given by the formula:
    # max_rank = min(d + r, n, d * r)
    
    # Calculate each term in the min function
    term1 = d + r
    term2 = n
    term3 = d * r
    
    # Calculate the final result
    max_rank = min(term1, term2, term3)
    
    # Print the equation and the result
    print("The maximal rank is determined by the formula: min(d + r, n, d * r)")
    print(f"Substituting the values:")
    print(f"max_rank = min({d} + {r}, {n}, {d} * {r})")
    print(f"max_rank = min({term1}, {term2}, {term3})")
    print(f"The maximal rank of the complementary Choi matrix is: {max_rank}")
    
    return max_rank

if __name__ == '__main__':
    # Example values for d, n, and r
    # Feel free to change these values to test different scenarios
    d_val = 5  # Dimension of H1
    n_val = 20 # Dimension of H2
    r_val = 3  # Rank of J(Lambda)
    
    calculate_max_complementary_rank(d_val, n_val, r_val)
    print("\n" + "="*40 + "\n")
    
    # Another example from the thought process
    d_val_2 = 2
    n_val_2 = 10
    r_val_2 = 1
    calculate_max_complementary_rank(d_val_2, n_val_2, r_val_2)
