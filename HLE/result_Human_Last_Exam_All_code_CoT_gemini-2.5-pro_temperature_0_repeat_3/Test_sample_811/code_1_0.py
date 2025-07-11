import math
from scipy.special import i0

def solve_probability():
    """
    Calculates the probability that the given matrix is diagonalisable.

    The probability is given by the formula:
    P(diag) = 1 - P(non-diag)
    P(non-diag) = e^(-1) * (sum_{k=1 to inf} 1/(k!)^3) / (sum_{k=0 to inf} 1/(k!)^2)
    """

    # The denominator sum, sum_{k=0 to inf} 1/(k!)^2, is I_0(2),
    # where I_0 is the modified Bessel function of the first kind of order 0.
    sum_denom = i0(2)

    # The numerator sum, sum_{k=1 to inf} 1/(k!)^3, is computed numerically.
    # We first compute sum_{k=0 to inf} 1/(k!)^3 and then subtract the k=0 term (which is 1).
    sum_num_from_0 = 0.0
    k = 0
    while True:
        try:
            term = 1.0 / (math.factorial(k)**3)
            if term < 1e-20:  # Stop when terms become negligible
                break
            sum_num_from_0 += term
            k += 1
        except OverflowError: # factorial grows very fast
            break
            
    sum_num = sum_num_from_0 - 1.0

    # The factor e^-1
    exp_minus_1 = math.exp(-1)

    # Probability of being non-diagonalisable
    prob_non_diag = exp_minus_1 * sum_num / sum_denom

    # Probability of being diagonalisable
    prob_diag = 1 - prob_non_diag

    print("The probability of the matrix being diagonalisable is calculated using the formula:")
    print("P(diag) = 1 - C1 * C2 / C3\n")
    print("The values of the components are:")
    print(f"C1 = e^-1 = {exp_minus_1}")
    print(f"C2 = sum_{k=1 to inf} 1/(k!)^3 = {sum_num}")
    print(f"C3 = sum_{k=0 to inf} 1/(k!)^2 = {sum_denom}\n")
    print(f"The probability of the matrix being non-diagonalisable is: {prob_non_diag}")
    print(f"The final probability of the matrix being diagonalisable is: {prob_diag}")
    
    return prob_diag

final_answer = solve_probability()
# The final answer is requested in a specific format.
# print(f"\n<<<{final_answer}>>>")