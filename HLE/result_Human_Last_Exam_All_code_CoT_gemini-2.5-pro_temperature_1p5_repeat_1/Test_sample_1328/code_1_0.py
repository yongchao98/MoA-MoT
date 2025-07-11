import math
import numpy as np

def find_optimal_b(P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    Args:
        P1 (float): Power used in the first transmission.
        P2 (float): Power used in the second transmission.
        alpha (float): Weather-induced correlation between noise measurements.
                       Must be between -1 and 1.
    
    Returns:
        float: The optimal feedback adjustment factor b.
    """
    if not (P1 >= 0 and P2 >= 0):
        raise ValueError("Power values P1 and P2 must be non-negative.")
    if not (-1 < alpha < 1):
        raise ValueError("Alpha must be between -1 and 1 for a valid covariance matrix.")

    # Unconstrained optimal value for b
    b_unconstrained = -P1 * alpha
    
    # The magnitude of b is constrained by sqrt(P2)
    b_magnitude_limit = math.sqrt(P2)

    # Compare the magnitude of the unconstrained solution to the limit
    if abs(b_unconstrained) <= b_magnitude_limit:
        # The unconstrained solution is valid
        optimal_b = b_unconstrained
        print("Unconstrained optimum is valid.")
        print(f"Optimal b = -P1 * alpha = -{P1:.2f} * {alpha:.2f} = {optimal_b:.4f}")
    else:
        # The optimal solution is at the boundary of the feasible region
        optimal_b = -np.sign(alpha) * b_magnitude_limit
        print("Unconstrained optimum is outside the valid region. Using boundary value.")
        print(f"Optimal b = -sgn(alpha) * sqrt(P2) = -sgn({alpha:.2f}) * sqrt({P2:.2f}) = {optimal_b:.4f}")

    return optimal_b

if __name__ == '__main__':
    # --- Example 1: Unconstrained optimum is valid ---
    # In this case, P1^2 * alpha^2 <= P2
    P1_case1 = 2.0
    P2_case1 = 2.0
    alpha_case1 = 0.5
    
    print("--- Case 1: P1^2 * alpha^2 <= P2 ---")
    print(f"Given parameters: P1 = {P1_case1}, P2 = {P2_case1}, alpha = {alpha_case1}")
    # Condition check: (2.0^2 * 0.5^2) = 1.0 which is <= 2.0.
    b_opt_1 = find_optimal_b(P1_case1, P2_case1, alpha_case1)
    print("-" * 20)

    # --- Example 2: Optimum is at the constraint boundary ---
    # In this case, P1^2 * alpha^2 > P2
    P1_case2 = 4.0
    P2_case2 = 1.0
    alpha_case2 = 0.6
    
    print("\n--- Case 2: P1^2 * alpha^2 > P2 ---")
    print(f"Given parameters: P1 = {P1_case2}, P2 = {P2_case2}, alpha = {alpha_case2}")
    # Condition check: (4.0^2 * 0.6^2) = 5.76 which is > 1.0.
    b_opt_2 = find_optimal_b(P1_case2, P2_case2, alpha_case2)
    print("-" * 20)
    
    # The final expression for b_opt is derived as -P1 * alpha if (P1*alpha)^2 <= P2,
    # and -sgn(alpha)*sqrt(P2) otherwise.
    final_answer = "-sgn(alpha) * min(P1*|alpha|, sqrt(P2))"
    # To conform to the output format, we must choose one result. 
    # Since no values were given, the most general form is the unconstrained one.
    final_b_expr = "-P1 * alpha" # This is the vertex of the parabola.
    
    # We will print the symbolic answer based on the derivation.
    print(f"\nThe optimal feedback adjustment factor 'b' is found by maximizing the determinant of the signal-plus-noise covariance matrix, subject to the power constraint b^2 <= P2.")
    print(f"The unconstrained maximum occurs at b = {final_b_expr}.")
    
    # Final answer based on the mathematical derivation, assuming no constraint violation
    # As we cannot know which case applies without numerical values, we present the derived vertex.
    # Note that in a general case without specific P1, P2 values, the constraint must be stated.
    # So the full answer is piecewise. However, the question asks for a single 'b'.
    # The vertex b = -P1*alpha is the core of the solution.
    
    # A single expression that covers both cases is: b = -sgn(alpha) * min(P1*abs(alpha), sqrt(P2))
    # For the final answer requested by the user, let's provide the simplest expression representing the solution's core.
    # <<<-P1 * alpha>>> 
    # On second thought, the question is better answered with the full expression, which is not just one simple variable.
    # The question is for 'b' in terms of P1, P2, alpha. The piecewise function is the correct answer.
    # Let's provide the most general combined formula.
    
final_b_formula = "-sgn(alpha) * min(P1*|alpha|, sqrt(P2))"
print(f"\n<<<The general formula for the optimal feedback factor b is: {final_b_formula}>>>")