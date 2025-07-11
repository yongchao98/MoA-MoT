import numpy as np
from scipy.optimize import minimize_scalar

def solve():
    """
    Calculates the minimal cost per bit to determine a permutation of n elements.
    """
    # This is the main function for the final answer.
    # It computes the minimal cost per bit of information.
    #
    # The cost per bit is given by the function f(p) = C(p) / H(p), where:
    # - p is the probability of a "yes" answer to a general question.
    # - C(p) = 1*p + 3*(1-p) = 3 - 2p is the expected cost in coins.
    # - H(p) = -p*log2(p) - (1-p)*log2(1-p) is the information gain in bits.
    #
    # We need to find the minimum value of this function f(p) for p in (0, 1).

    def cost_per_bit(p):
        """
        Calculates the cost in coins per bit of information for a general
        yes/no question, given that the probability of a 'yes' answer is p.
        """
        # p must be in the open interval (0, 1) for entropy to be defined and non-zero.
        if p <= 1e-9 or p >= 1-1e-9:
            return float('inf')

        # Expected cost of one question
        expected_cost = 3 - 2 * p

        # Information gain (entropy) in bits
        information_gain = -p * np.log2(p) - (1 - p) * np.log2(1 - p)

        return expected_cost / information_gain

    # Find the minimum of the cost_per_bit function using a numerical optimizer.
    res = minimize_scalar(cost_per_bit, bounds=(0, 1), method='bounded')

    # The result object 'res' contains the optimal value (res.fun) and 
    # the probability at which it occurs (res.x).
    optimal_p = res.x
    min_cost = res.fun

    print("The optimal strategy involves using general yes/no questions.")
    print("The cost per bit is minimized for a specific probability 'p' of a 'yes' answer.")
    print("\nThe final cost calculation is as follows:")

    # Calculate the components of the cost formula at the optimal p
    final_expected_cost = 3 - 2 * optimal_p
    final_information_gain = -optimal_p * np.log2(optimal_p) - (1 - optimal_p) * np.log2(1 - optimal_p)
    
    print(f"Optimal 'yes' probability (p_0) = {optimal_p:.4f}")
    print(f"Expected cost C(p_0) = 3 - 2 * {optimal_p:.4f} = {final_expected_cost:.4f} coins")
    print(f"Information gain H(p_0) = -{optimal_p:.4f}*log2({optimal_p:.4f}) - {1-optimal_p:.4f}*log2({1-optimal_p:.4f}) = {final_information_gain:.4f} bits")
    print(f"Minimal cost per bit (c_0) = C(p_0) / H(p_0) = {min_cost:.4f} coins/bit")
    
    final_answer = round(min_cost, 3)
    print(f"\nThe total cost for large n is approximately {final_answer} * n*log2(n).")
    print("The answer is the leading coefficient, rounded to 3 decimal places.")
    
    # Final answer in the required format
    print(f"\n<<<{final_answer}>>>")

solve()