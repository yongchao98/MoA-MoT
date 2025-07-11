import math

def find_optimal_b(P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor b that maximizes mutual information.

    The optimal value for b is derived by maximizing the determinant of the
    output covariance matrix, |K_{X+Z}|, subject to power constraints.
    The function to maximize with respect to b is f(b) = -b^2 - 2*alpha*P1*b.
    This is a downward-opening parabola with a vertex at b = -alpha * P1.

    The power constraint on the second transmission implies b^2 <= P2.

    Args:
        P1 (float): Power used in the first transmission. Must be non-negative.
        P2 (float): Power used in the second transmission. Must be non-negative.
        alpha (float): Noise correlation coefficient. Must be in [-1, 1].

    Returns:
        float: The optimal feedback factor b.
        str: A string explaining the calculation.
    """
    # Input validation
    if P1 < 0 or P2 < 0:
        raise ValueError("Power values P1 and P2 cannot be negative.")
    # For K to be a valid covariance matrix, |alpha| must be <= 1.
    if not -1 <= alpha <= 1:
        raise ValueError("Correlation alpha must be between -1 and 1.")

    # The unconstrained optimal value for b is at the vertex of the parabola.
    b_unconstrained = -alpha * P1

    # Check if the unconstrained optimum is within the feasible region defined by b^2 <= P2.
    if b_unconstrained**2 <= P2:
        # The ideal value is feasible.
        optimal_b = b_unconstrained
        calculation_str = f"b = -({alpha}) * {P1} = {optimal_b:.4f}"
    else:
        # The ideal value is not feasible, so we saturate at the boundary of the
        # feasible region [-sqrt(P2), sqrt(P2)] that is closest to the vertex.
        # This corresponds to choosing the sign that is opposite to alpha's sign.
        if alpha == 0:
            optimal_b = 0.0
            calculation_str = f"b = 0.0, since alpha is 0"
        else:
            optimal_b = -math.copysign(math.sqrt(P2), alpha)
            calculation_str = f"b = -sign({alpha}) * sqrt({P2}) = {optimal_b:.4f}"
            
    return optimal_b, calculation_str

# --- Example Usage ---
# You can change these values to test different scenarios.

# Scenario 1: The ideal feedback value is within the power limits.
print("--- Scenario 1: Unconstrained Case ---")
P1_case1, P2_case1, alpha_case1 = 3.0, 10.0, 0.5
print(f"Given parameters: P1 = {P1_case1}, P2 = {P2_case1}, alpha = {alpha_case1}")
b_ideal_sq_1 = (alpha_case1 * P1_case1)**2
print(f"Condition check: (alpha*P1)^2 = {b_ideal_sq_1:.4f}. P2 = {P2_case1}.")
print(f"Since {b_ideal_sq_1:.4f} <= {P2_case1}, the ideal value is feasible.")
optimal_b1, calc1_str = find_optimal_b(P1_case1, P2_case1, alpha_case1)
print(f"Optimal b Calculation: {calc1_str}")

# Scenario 2: The ideal feedback requires more power than available.
print("\n--- Scenario 2: Constrained Case ---")
P1_case2, P2_case2, alpha_case2 = 8.0, 10.0, 0.5
print(f"Given parameters: P1 = {P1_case2}, P2 = {P2_case2}, alpha = {alpha_case2}")
b_ideal_sq_2 = (alpha_case2 * P1_case2)**2
print(f"Condition check: (alpha*P1)^2 = {b_ideal_sq_2:.4f}. P2 = {P2_case2}.")
print(f"Since {b_ideal_sq_2:.4f} > {P2_case2}, the solution is constrained by P2.")
optimal_b2, calc2_str = find_optimal_b(P1_case2, P2_case2, alpha_case2)
print(f"Optimal b Calculation: {calc2_str}")