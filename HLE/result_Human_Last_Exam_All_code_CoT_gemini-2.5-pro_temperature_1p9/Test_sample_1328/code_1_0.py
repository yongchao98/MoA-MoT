import math

def find_optimal_b(P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    Args:
        P1 (float): Power used in the first transmission.
        P2 (float): Power used in the second transmission.
        alpha (float): Weather-induced correlation between noise measurements.

    Returns:
        float: The optimal value of b.
    """
    print("--- Calculation of Optimal Feedback Factor b ---")
    
    # Validate inputs
    if P1 < 0 or P2 < 0:
        raise ValueError("Power values P1 and P2 cannot be negative.")
    if not -1 < alpha < 1:
        raise ValueError("Correlation alpha must be between -1 and 1.")

    print("\nGiven parameters:")
    print(f"P1 = {P1}")
    print(f"P2 = {P2}")
    print(f"alpha = {alpha}")

    # Step 1: Calculate the unconstrained optimal b
    b_unc = -alpha * P1
    print("\n1. Calculate the unconstrained optimal b:")
    print(f"   b_unc = -alpha * P1 = -{alpha} * {P1} = {b_unc}")

    # Step 2: Calculate the feasible range for b from the power constraint P_U2 >= 0
    # b^2 <= P2  => |b| <= sqrt(P2)
    b_boundary = math.sqrt(P2)
    print("\n2. Calculate the feasible range for b, |b| <= sqrt(P2):")
    print(f"   The boundary value is sqrt(P2) = sqrt({P2}) = {b_boundary:.4f}")
    print(f"   The feasible interval for b is [{-b_boundary:.4f}, {b_boundary:.4f}]")

    # Step 3: Apply the constraint to find the final optimal b
    # This is equivalent to clipping b_unc to the interval [-b_boundary, b_boundary]
    if b_unc > b_boundary:
        b_optimal = b_boundary
        clipping_reason = f"b_unc ({b_unc:.4f}) > boundary ({b_boundary:.4f}), so b is clipped to the upper bound."
    elif b_unc < -b_boundary:
        b_optimal = -b_boundary
        clipping_reason = f"b_unc ({b_unc:.4f}) < -boundary (-{b_boundary:.4f}), so b is clipped to the lower bound."
    else:
        b_optimal = b_unc
        clipping_reason = f"b_unc ({b_unc:.4f}) is within the feasible range, no clipping needed."
        
    print("\n3. Apply the constraint to find the final optimal b:")
    print(f"   {clipping_reason}")
    print("\n--- Final Result ---")
    print(f"The optimal feedback adjustment factor b is: {b_optimal}")
    
    return b_optimal

# --- Example Usage ---
# You can change these values to see how the result changes.

# Example 1: Unconstrained optimum is inside the feasible range
print("---------- Example 1 ----------")
find_optimal_b(P1=4, P2=10, alpha=0.5)

# Example 2: Unconstrained optimum is outside the feasible range
print("\n\n---------- Example 2 ----------")
find_optimal_b(P1=8, P2=9, alpha=0.6)