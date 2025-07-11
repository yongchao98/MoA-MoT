import math

def solve_diagonalisable_probability():
    """
    Calculates the probability that the random matrix is diagonalisable.
    The method involves finding the probability of the complementary event (non-diagonalisability)
    using a recursive formulation based on the i.i.d. nature of the Poisson variables.
    This leads to a formula for the probability of non-diagonalisability, Q, as a ratio of two infinite sums.
    The script numerically computes these sums to find the final answer, P = 1 - Q.
    """

    # The sums are approximated by truncating the series at K.
    # The terms p_k = exp(-1)/k! decrease very rapidly, so K=30 is highly accurate.
    K = 30

    # Pre-calculate probabilities p_k = P(X=k) for a Poisson(1) distribution.
    p = [math.exp(-1) / math.factorial(k) for k in range(K + 1)]

    # Calculate the numerator sum D3 from the derived formula for Q.
    # D3 = sum_{k=1 to inf} p_k^3 / (1 + p_k)
    numerator_sum = sum((p[k]**3) / (1 + p[k]) for k in range(1, K + 1))

    # Calculate the denominator sum D2 from the derived formula for Q.
    # D2 = sum_{k=0 to inf} p_k^2 / (1 + p_k)
    denominator_sum = sum((p[k]**2) / (1 + p[k]) for k in range(K + 1))

    # Q is the probability of the matrix being non-diagonalisable.
    prob_non_diag = numerator_sum / denominator_sum

    # The final answer is the probability of the matrix being diagonalisable.
    prob_diag = 1 - prob_non_diag

    print("The probability of the matrix being diagonalisable is P = 1 - Q,")
    print("where Q is the probability of it being non-diagonalisable.")
    print("The formula for Q is the ratio of two sums, D3 / D2.")
    print("\nNumerically evaluating the components of the equation:")
    print(f"Numerator Sum (D3) = {numerator_sum}")
    print(f"Denominator Sum (D2) = {denominator_sum}")
    print("\nThe equation for the probability of being non-diagonalisable is:")
    print(f"Q = {numerator_sum} / {denominator_sum} = {prob_non_diag}")
    print("\nFinally, the equation for the probability of being diagonalisable is:")
    print(f"P = 1 - {prob_non_diag} = {prob_diag}")
    
    # Return the final probability for the answer tag.
    return prob_diag

# Execute the calculation and store the result.
final_probability = solve_diagonalisable_probability()
print(f"\n<<<>>>\n{final_probability}\n<<<>>>")