import math

def calculate_optimal_b(p1, p2, alpha):
    """
    Calculates the optimal feedback adjustment factor b.

    Args:
        p1 (float): Power used in the first transmission.
        p2 (float): Power used in the second transmission.
        alpha (float): Weather-induced correlation.
    """
    # The optimal feedback factor b is found by maximizing the determinant of K_{X+Z}
    # The derivative of the determinant with respect to b set to zero gives:
    # b = -p1 * alpha
    b_optimal = -p1 * alpha

    # We must also respect the physical power constraint where b^2 <= p2.
    # The following example parameters are chosen to satisfy this.
    # We check the condition:
    is_valid = (b_optimal**2) <= p2

    print(f"Given Parameters:")
    print(f"  Power P1 = {p1}")
    print(f"  Power P2 = {p2}")
    print(f"  Correlation alpha = {alpha}\n")

    print("Step 1: The optimal feedback factor is derived from maximizing the channel capacity.")
    print("The result of this maximization is the formula: b = -P1 * alpha.\n")

    print("Step 2: Substitute the given values into the formula:")
    print(f"b = -{p1} * {alpha}")
    b_final = -p1 * alpha
    print(f"b = {b_final:.4f}\n")
    
    print("Step 3: Check the physical power constraint (b^2 <= P2).")
    print(f"  b^2 = ({b_final:.4f})^2 = {b_final**2:.4f}")
    print(f"  Is {b_final**2:.4f} <= {p2}? {'Yes' if is_valid else 'No'}")
    if is_valid:
        print("The constraint is satisfied, so this is the final optimal value for b.")
    else:
        print("Warning: Constraint is not satisfied. The optimal b would be on the boundary.")


# --- Example Calculation ---
# Let's use example values for P1, P2, and alpha.
# These values are selected such that the physical constraint is met.
P1_val = 5.0
P2_val = 10.0
alpha_val = 0.6

calculate_optimal_b(P1_val, P2_val, alpha_val)
