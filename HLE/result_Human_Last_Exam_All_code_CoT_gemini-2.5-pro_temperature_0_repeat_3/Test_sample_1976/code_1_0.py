import math

def calculate_norm_t1(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n.

    This function implements the full summation formula derived from the
    Bloch decomposition and also prints the simplified analytical result.

    Args:
        n: An odd integer.
    """
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print("Error: n must be a positive odd integer.")
        return

    N = n + 1

    # Numerator is the sum: sum_{m=1}^{N} C(N, m) * 3^m * |1 + (-1)^m * 3^(n-m)|
    numerator_sum = 0
    for m in range(1, N + 1):
        # Use math.comb for combinations C(N, m)
        comb = math.comb(N, m)
        
        # Calculate the term inside the absolute value
        inner_term = 1 + ((-1)**m) * (3**(n - m))
        
        # Add the contribution of the m-th term to the sum
        term = comb * (3**m) * abs(inner_term)
        numerator_sum += term

    # Denominator is (1 + 3^n) * 2^N
    denominator = (1 + 3**n) * (2**N)

    # The 1-norm is the ratio of the numerator and denominator
    norm_t1 = numerator_sum / denominator

    # The final equation for the 1-norm
    print(f"For n = {n}:")
    print(f"The numerator of the expression is: {numerator_sum}")
    print(f"The denominator of the expression is: {denominator}")
    print(f"The 1-norm of the correlation matrix T is: {norm_t1}")

    # Verification using the simplified formula: 1 - 1/2^(n+1)
    simplified_result = 1 - 1 / (2**(n + 1))
    print(f"\nVerification with the simplified formula (1 - 1/2^(n+1)): {simplified_result}")
    # Check if the results are close enough to account for potential floating point inaccuracies
    if math.isclose(norm_t1, simplified_result):
        print("The calculated result matches the simplified formula.")
    else:
        print("The calculated result does not match the simplified formula.")


if __name__ == '__main__':
    # You can change this value to any positive odd integer
    odd_n = 3
    calculate_norm_t1(odd_n)