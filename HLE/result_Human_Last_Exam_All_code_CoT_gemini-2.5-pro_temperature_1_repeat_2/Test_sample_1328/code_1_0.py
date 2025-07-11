import math

def find_optimal_b(alpha, p1, p2):
    """
    Calculates the optimal feedback adjustment factor 'b' to maximize mutual information.

    Args:
        alpha (float): The weather-induced correlation between noise measurements.
        p1 (float): Power used in the first transmission.
        p2 (float): Power used in the second transmission.

    Returns:
        float: The optimal value of the feedback adjustment factor 'b'.
    """
    if p1 < 0 or p2 < 0:
        print("Error: Power values (P1, P2) must be non-negative.")
        return None
    
    # Step 1: Find the unconstrained optimal 'b' by maximizing det(K_Y).
    # The determinant is a downward-opening parabola in 'b': det(K_Y) = C - b^2 - 2*b*alpha*P1
    # The maximum is at the vertex, where the derivative is zero:
    # d/db(det(K_Y)) = -2*b - 2*alpha*P1 = 0  => b = -alpha * P1
    b_unconstrained = -alpha * p1

    # Step 2: Consider the physical constraint on 'b'.
    # The power of the information-bearing part of the second signal (P_X_tilde_2) must be non-negative.
    # P2 = E[X2^2] = E[(X_tilde_2 - b*Z1)^2] = P_X_tilde_2 + b^2.
    # So, P_X_tilde_2 = P2 - b^2 >= 0, which implies b^2 <= P2.
    # This means 'b' must be in the interval [-sqrt(P2), sqrt(P2)].
    b_limit = math.sqrt(p2)

    # Step 3: Find the optimal 'b' within the allowed interval.
    # The optimal value is the point in the interval [-b_limit, b_limit]
    # that is closest to the unconstrained optimum.
    
    print(f"Given parameters: alpha = {alpha}, P1 = {p1}, P2 = {p2}")
    print(f"The unconstrained optimal value for b is b = -alpha * P1 = -({alpha}) * ({p1}) = {b_unconstrained:.4f}")
    print(f"The constraint from P2 is |b| <= sqrt(P2), so |b| <= {b_limit:.4f}")

    if abs(b_unconstrained) <= b_limit:
        optimal_b = b_unconstrained
        print("\nThe unconstrained optimum is within the allowed range.")
        print(f"Final optimal b = {optimal_b:.4f}")
    else:
        # If the unconstrained optimum is outside the interval, the optimum is the boundary point
        # closest to it. This can be expressed by clipping the value.
        optimal_b = math.copysign(b_limit, b_unconstrained)
        print("\nThe unconstrained optimum is outside the allowed range.")
        print(f"The optimal b is clipped to the boundary.")
        print(f"b = sign({b_unconstrained:.4f}) * min(|{b_unconstrained:.4f}|, {b_limit:.4f})")
        print(f"b = {math.copysign(1, b_unconstrained):.0f} * {b_limit:.4f} = {optimal_b:.4f}")

    return optimal_b

# --- Example Usage ---
# Case 1: Unconstrained optimum is feasible
print("--- Case 1: Unconstrained optimum is feasible ---")
find_optimal_b(alpha=0.5, p1=5, p2=10)
print("\n" + "="*50 + "\n")

# Case 2: Unconstrained optimum is not feasible
print("--- Case 2: Unconstrained optimum is not feasible ---")
find_optimal_b(alpha=0.5, p1=10, p2=12)
