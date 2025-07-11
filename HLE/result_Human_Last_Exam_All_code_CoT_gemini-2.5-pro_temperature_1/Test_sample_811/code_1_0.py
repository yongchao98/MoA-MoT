import math

def calculate_diagonalizable_probability():
    """
    Calculates the probability that the matrix is diagonalisable.
    """
    # The terms in the sums decrease very rapidly due to the k! term.
    # A limit of K=30 is more than sufficient for high precision.
    K = 30
    
    # Calculate p_k = exp(-1) / k! for k = 0, 1, ..., K
    p = [0.0] * (K + 1)
    e_inv = math.exp(-1)
    for k in range(K + 1):
        p[k] = e_inv / math.factorial(k)

    # Calculate the denominator: S = sum_{j=0 to inf} p_j^2 / (1 + p_j)
    denominator_sum = 0
    for j in range(K + 1):
        denominator_sum += (p[j]**2) / (1 + p[j])

    # Calculate the numerator of the non-diagonalisable part:
    # N = sum_{k=1 to inf} p_k^3 / (1 + p_k)
    numerator_sum = 0
    for k in range(1, K + 1):
        numerator_sum += (p[k]**3) / (1 + p[k])

    # The probability of being non-diagonalisable is Numerator / Denominator
    prob_non_diagonalisable = numerator_sum / denominator_sum
    
    # The probability of being diagonalisable is 1 - prob_non_diagonalisable
    prob_diagonalisable = 1 - prob_non_diagonalisable

    # Print the final equation with the computed values
    print("The probability is calculated using the formula:")
    print("P(diagonalisable) = 1 - [sum_{k=1 to inf} p_k^3 / (1+p_k)] / [sum_{j=0 to inf} p_j^2 / (1+p_j)]")
    print("\nWith numerical evaluation:")
    print(f"Numerator sum = {numerator_sum}")
    print(f"Denominator sum = {denominator_sum}")
    print(f"P(diagonalisable) = 1 - {numerator_sum} / {denominator_sum}")
    print(f"P(diagonalisable) = {prob_diagonalisable}")

if __name__ == '__main__':
    calculate_diagonalizable_probability()
