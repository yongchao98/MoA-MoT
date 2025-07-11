import math
import scipy.optimize

def solve_and_calculate_cost_factor():
    """
    Solves for the optimal cost factor for the sorting problem.

    The problem reduces to finding a constant k such that the total cost is k * log2(n!).
    This constant k is the minimum cost per bit of information. We derived that k can be
    found by solving the equation p^3 + p - 1 = 0 and then calculating k = 1 / (-log2(p)).
    """

    # The equation to solve for the optimal probability split 'p' is p^3 + p - 1 = 0.
    # The coefficients of the polynomial p^3 + 0*p^2 + 1*p - 1 are 1, 0, 1, -1.
    c_p3 = 1
    c_p1 = 1
    c_const = -1
    
    print(f"The analysis leads to the equation: {c_p3}*p^3 + {c_p1}*p + {c_const} = 0")

    # Define the function for the root-finding algorithm.
    def poly(p):
        return p**3 + p - 1

    # Find the real root of the polynomial in the interval (0, 1).
    # A root must exist here since poly(0) = -1 and poly(1) = 1.
    try:
        p_root = scipy.optimize.brentq(poly, 0, 1)
        print(f"The optimal probability split p₀ is approximately: {p_root:.6f}")
    except (ImportError, ModuleNotFoundError):
        print("\nSciPy is not installed. Using a simplified numerical approach (Newton's method).")
        # Fallback to a simple Newton's method if scipy is not available
        p_root = 0.6  # Initial guess
        for _ in range(5):
             p_root = p_root - poly(p_root) / (3 * p_root**2 + 1)
        print(f"The optimal probability split p₀ is approximately: {p_root:.6f}")


    # The minimal cost per bit of information, k, is 1 / (-log2(p_root)).
    cost_factor = 1 / (-math.log2(p_root))
    
    print(f"The minimal cost per bit of information is: {cost_factor:.6f}")
    print("\nThe final answer is this cost factor, rounded to 3 decimal places.")
    
    # Final formatted output
    final_answer = f"{cost_factor:.3f}"
    
    return final_answer

final_answer_value = solve_and_calculate_cost_factor()
print(f"<<<{final_answer_value}>>>")
