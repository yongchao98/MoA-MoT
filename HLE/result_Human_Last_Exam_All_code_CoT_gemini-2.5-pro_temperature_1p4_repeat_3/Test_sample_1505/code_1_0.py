import math

def approximate_sum(n: int):
    """
    Calculates the approximation for the sum S(n) for a given positive integer n.
    The formula used is n^2/2 + 1/120 + 1/(252*n), which has an
    absolute error of O(n^-2).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    # The derived approximation formula is n^2/2 + 1/120 + 1/(252*n).
    # Let's calculate the value of each term in this formula.
    term1_val = n**2 / 2
    term2_val = 1 / 120
    term3_val = 1 / (252 * n)

    # The result is the sum of these terms.
    approximation_result = term1_val + term2_val + term3_val

    # As requested, we print the final equation with each number.
    # We format the floating point numbers for better readability.
    print(f"For n = {n}, the sum is approximated by the formula: n^2/2 + 1/120 + 1/(252*n)")
    print("The final equation with the calculated values for each term is:")
    print(f"{term1_val} + {term2_val:.8f} + {term3_val:.8f} = {approximation_result:.8f}")


# Example usage with n = 100. The user can change this value.
n_value = 100
approximate_sum(n_value)
