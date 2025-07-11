import math

def find_optimal_feedback_factor():
    """
    Calculates the optimal feedback adjustment factor 'b' for a wireless channel
    with memory to maximize mutual information.
    """

    # --- Problem Parameters ---
    # You can change these values to reflect your specific system.
    # Power used in the first transmission
    P1 = 5.0
    # Power used in the second transmission
    P2 = 4.0
    # Weather-induced correlation between noise measurements
    alpha = 0.6

    print(f"Given parameters:")
    print(f"Power in first transmission (P1): {P1}")
    print(f"Power in second transmission (P2): {P2}")
    print(f"Noise correlation (alpha): {alpha}\n")

    # --- Step 1: Derive the unconstrained optimal 'b' ---
    # The mutual information is maximized by maximizing the determinant of the 
    # output covariance matrix, |K_{X+Z}|.
    # This determinant can be expressed as a quadratic function of b:
    # |K_{X+Z}| = -b^2 + (2 * P1 * alpha) * b + Constant
    #
    # To find the maximum, we take the derivative with respect to b and set it to 0:
    # d/db |K_{X+Z}| = -2*b + 2*P1*alpha = 0
    # This gives the unconstrained solution for b.
    b_unconstrained = P1 * alpha

    print("--- Calculation Steps ---")
    print(f"1. The unconstrained optimal 'b' that maximizes the objective function is b = P1 * alpha.")
    print(f"   b_unconstrained = {P1} * {alpha} = {b_unconstrained}")

    # --- Step 2: Apply the physical power constraint ---
    # The total power in the second slot is P2 = E[X2^2].
    # Since X2 = S2 + b*Z1 (where S2 is the new signal), the power of the new signal is
    # P_S2 = P2 - b^2. Physical signals cannot have negative power, so P_S2 >= 0.
    # This implies the constraint: b^2 <= P2.
    
    print("\n2. Applying the physical constraint b^2 <= P2:")
    print(f"   ({b_unconstrained})^2 <= {P2}  =>  {b_unconstrained**2} <= {P2}?")

    # --- Step 3: Determine the final optimal 'b' ---
    # Check if the unconstrained solution is physically feasible.
    if b_unconstrained**2 <= P2:
        # The unconstrained optimum is feasible.
        b_optimal = b_unconstrained
        print("   The condition is satisfied. The optimal 'b' is the unconstrained solution.")
        print("\nFinal Equation:")
        print(f"b_optimal = P1 * alpha")
        print(f"b_optimal = {P1} * {alpha} = {b_optimal}")

    else:
        # The unconstrained optimum is not feasible. The optimal value must lie
        # on the boundary of the feasible region, [-sqrt(P2), sqrt(P2)].
        # We choose the boundary point that is closest to the unconstrained optimum.
        b_optimal = math.copysign(math.sqrt(P2), b_unconstrained)
        print("   The condition is NOT satisfied.")
        print("   The optimal 'b' is clipped to the boundary of the feasible region [-sqrt(P2), sqrt(P2)].")
        print("\nFinal Equation:")
        print(f"b_optimal = sign(P1 * alpha) * sqrt(P2)")
        print(f"b_optimal = sign({P1} * {alpha}) * sqrt({P2}) = {b_optimal}")

    print("\n-------------------------------------------")
    print(f"The final optimal feedback adjustment factor b is: {b_optimal}")
    print("-------------------------------------------")

    return b_optimal

if __name__ == '__main__':
    optimal_b = find_optimal_feedback_factor()
    # The final answer format as requested.
    # print(f"<<<{optimal_b}>>>") 
    
find_optimal_feedback_factor()
<<<2.0>>>