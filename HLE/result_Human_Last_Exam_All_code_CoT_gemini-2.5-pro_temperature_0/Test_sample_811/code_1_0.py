import math

def solve_probability():
    """
    Calculates the probability that the given matrix is diagonalisable.
    """
    # The matrix is not diagonalisable if and only if X_N = X_{N+2} and X_N != 0.
    # Let pi_ND be the probability of this event. The answer is 1 - pi_ND.
    # The probability pi_ND can be expressed as a ratio of two infinite sums involving
    # p_k = P(Poisson(1) = k) = exp(-1) / k!
    #
    # pi_ND = (sum_{k=1 to inf} p_k^3 / (1+p_k)) / (sum_{k=0 to inf} p_k^2 / (1+p_k))
    #
    # We will numerically calculate this value by truncating the infinite sums.
    # The terms decrease very rapidly, so a small number of terms is sufficient.

    K_MAX = 30  # Truncation limit for the sums.

    # Pre-calculate log factorials for efficiency
    log_factorials = [0.0] * (K_MAX + 1)
    for i in range(1, K_MAX + 1):
        log_factorials[i] = log_factorials[i-1] + math.log(i)

    def get_p(k):
        """Computes p_k = e^(-1) / k! using log factorials for numerical stability."""
        if not (0 <= k < len(log_factorials)):
            return 0
        log_pk = -1 - log_factorials[k]
        return math.exp(log_pk)

    # Calculate the numerator sum for the probability of non-diagonalisability
    numerator_sum = 0
    for k in range(1, K_MAX):
        pk = get_p(k)
        if pk == 0:  # Further terms will be zero
            break
        numerator_sum += (pk**3) / (1 + pk)

    # Calculate the denominator sum
    denominator_sum = 0
    for k in range(K_MAX):
        pk = get_p(k)
        if pk == 0 and k > 1: # Further terms will be zero
            break
        denominator_sum += (pk**2) / (1 + pk)

    # Calculate the probability of being non-diagonalisable
    pi_ND = numerator_sum / denominator_sum
    
    # The final answer is the probability of being diagonalisable
    prob_diag = 1 - pi_ND

    print("The probability of the matrix being diagonalisable is P(diag) = 1 - P(non-diag).")
    print("P(non-diag) is calculated as the ratio of two sums derived from a Markov chain model.")
    print("\nLet p_k = e^(-1) / k! be the Poisson probability.")
    print("\nNumerator = Sum_{k=1 to inf} [ p_k^3 / (1 + p_k) ]")
    print(f"The approximated value of the numerator sum is: {numerator_sum}")
    print("\nDenominator = Sum_{k=0 to inf} [ p_k^2 / (1 + p_k) ]")
    print(f"The approximated value of the denominator sum is: {denominator_sum}")
    print("\nThe probability of the matrix being non-diagonalisable is the ratio:")
    print(f"P(non-diag) = {numerator_sum} / {denominator_sum} = {pi_ND}")
    print("\nFinally, the probability of the matrix being diagonalisable is:")
    print(f"P(diag) = 1 - {pi_ND} = {prob_diag}")

solve_probability()