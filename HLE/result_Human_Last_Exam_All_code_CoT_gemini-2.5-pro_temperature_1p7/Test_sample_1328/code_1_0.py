import math

def find_optimal_b(P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    Args:
        P1 (float): Power used in the first transmission.
        P2 (float): Power used in the second transmission.
        alpha (float): Weather-induced correlation between noise measurements.
                       Should be between -1 and 1.

    Returns:
        float: The optimal value of b.
    """
    print(f"Given parameters:")
    print(f"P1 = {P1}")
    print(f"P2 = {P2}")
    print(f"alpha = {alpha}\n")
    
    if not -1 <= alpha <= 1:
        print("Error: alpha must be between -1 and 1.")
        return None
    if P1 < 0 or P2 < 0:
        print("Error: Power values P1 and P2 cannot be negative.")
        return None

    # This is the square of the unconstrained optimal b
    unconstrained_b_squared = (P1 * alpha)**2

    print("Checking the condition: (P1 * alpha)^2 <= P2")
    print(f"({P1:.2f} * {alpha:.2f})^2 = {unconstrained_b_squared:.4f}")
    print(f"P2 = {P2:.4f}")

    if unconstrained_b_squared <= P2:
        # Case 1: The unconstrained optimum is feasible.
        b_opt = -P1 * alpha
        print("\nCondition is met. The optimal b is calculated as -P1 * alpha.")
        print(f"b_opt = -({P1:.2f}) * {alpha:.2f}")
        print(f"b_opt = {b_opt:.4f}")
    else:
        # Case 2: The optimum is on the boundary of the feasible region b^2 <= P2.
        b_opt = -math.copysign(1, alpha) * math.sqrt(P2)
        print(f"\nCondition is NOT met. The optimal b is on the boundary.")
        print("b_opt is calculated as -sign(alpha) * sqrt(P2).")
        sign_alpha = math.copysign(1, alpha)
        sqrt_P2 = math.sqrt(P2)
        print(f"b_opt = -sign({alpha:.2f}) * sqrt({P2:.2f}) = -({sign_alpha:.0f}) * {sqrt_P2:.4f}")
        print(f"b_opt = {b_opt:.4f}")
    
    return b_opt

# --- Example Usage ---
# You can change these values to test different scenarios.
P1_val = 3.0
P2_val = 1.0
alpha_val = 0.5
# Note: P1+P2 should not exceed the total budget 2P, but P is not needed for the calculation of b.

find_optimal_b(P1_val, P2_val, alpha_val)
