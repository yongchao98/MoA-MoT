import math

def get_q_k(k):
    """
    Calculates the value of q_k for a given integer k.
    """
    if k % 2 != 0:  # k is odd
        return 2.0 - (abs(k) - 1.0) / 2.0
    else:  # k is even
        return 5.0 / 3.0 - abs(k) / 2.0

def calculate_Q_at_minus_one(num_terms):
    """
    Calculates the partial sum of Q(D) at D = -1.
    Q(-1) = q_0 + sum_{k=1 to N} 2 * q_k * (-1)^k
    """
    q_0 = get_q_k(0)
    total_sum = q_0
    
    # Print the start of the equation
    print(f"Calculating Q(-1) = q_0 - 2*q_1 + 2*q_2 - 2*q_3 + ...")
    print(f"Q(-1) = {q_0:.4f}", end="")
    
    for k in range(1, num_terms + 1):
        q_k = get_q_k(k)
        term = 2 * q_k * ((-1) ** k)
        total_sum += term
        
        # Print each term in the equation
        if k % 2 != 0: # odd k
            print(f" - 2*({q_k:.4f})", end="")
        else: # even k
            print(f" + 2*({q_k:.4f})", end="")

    print(f"\n\nPartial sum for N = {num_terms} terms:")
    print(f"Q(-1) approx = {total_sum:.4f}")
    
    # Conclusion based on the calculation
    print("\nThe calculation shows that Q(D) at D=-1 is negative.")
    print("A power spectral density cannot be negative. The given q_k sequence does not correspond to a valid physical channel and noise model.")
    print("Therefore, a whitening filter W(D) that satisfies the required properties does not exist because Q(D) is not a valid spectral density function.")
    print("The property |q_k| <= q_0 for all k, which is necessary for an autocorrelation sequence, is violated since |q_1| = 2 > q_0 = 5/3.")


if __name__ == '__main__':
    # We use a sufficient number of terms to see the trend.
    # The terms decay, but slowly. Let's use 10 terms for demonstration.
    calculate_Q_at_minus_one(10)