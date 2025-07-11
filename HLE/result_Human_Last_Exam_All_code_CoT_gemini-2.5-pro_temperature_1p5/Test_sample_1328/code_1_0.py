import math

def solve_optimal_b(P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor b.

    Args:
        P1 (float): Power used in the first transmission. E[X1^2]
        P2 (float): Power used in the second transmission. E[X2^2]
        alpha (float): Weather-induced correlation between noise measurements.

    Returns:
        float: The optimal feedback adjustment factor b.
    """
    # The objective is to maximize |K_Y|, which simplifies to maximizing the function
    # g(b) = -b^2 - 2*P1*alpha*b
    # The constraint on b is derived from the fact that the power of the information
    # signal U2 must be non-negative: E[U2^2] = P2 - b^2 >= 0, so b^2 <= P2.
    # This means b must be in the interval [-sqrt(P2), sqrt(P2)].

    # The function g(b) is a downward-opening parabola. Its vertex is at b_v = -P1 * alpha.
    
    # We must find the maximum of this parabola on the interval [-sqrt(P2), sqrt(P2)].
    # The maximum occurs at the vertex if the vertex is within the interval.
    # Otherwise, it occurs at the boundary point closest to the vertex.

    # Condition for the vertex to be in the interval: (b_v)^2 <= P2
    # (-P1 * alpha)^2 <= P2  =>  P1^2 * alpha^2 <= P2
    
    print(f"Given parameters: P1 = {P1}, P2 = {P2}, alpha = {alpha}")

    if P1**2 * alpha**2 <= P2:
        # The vertex is within the feasible interval, so it is the optimum.
        b_optimal = -P1 * alpha
        print(f"Condition P1^2 * alpha^2 <= P2 is met ({P1**2 * alpha**2:.2f} <= {P2}).")
        print("The optimal b is at the vertex of the parabola.")
        print(f"b = -P1 * alpha = -({P1}) * ({alpha})")
    else:
        # The vertex is outside the interval. The optimum is the boundary point
        # closest to the vertex b_v = -P1 * alpha.
        # This is -sqrt(P2) if b_v < -sqrt(P2) (i.e., alpha > 0),
        # and +sqrt(P2) if b_v > +sqrt(P2) (i.e., alpha < 0).
        # This can be written compactly as -sign(alpha) * sqrt(P2).
        
        # Note: If alpha is 0, sign(alpha) is 0, so b is 0.
        # If P2 is 0, b is 0.
        sign_alpha = 0
        if alpha > 0:
            sign_alpha = 1
        elif alpha < 0:
            sign_alpha = -1
            
        b_optimal = -sign_alpha * math.sqrt(P2)
        print(f"Condition P1^2 * alpha^2 <= P2 is not met ({P1**2 * alpha**2:.2f} > {P2}).")
        print("The optimal b is at the boundary of the feasible region.")
        if sign_alpha != 0:
            print(f"b = -sign(alpha) * sqrt(P2) = -({sign_alpha}) * sqrt({P2})")
        else:
            print(f"b = -sign(alpha) * sqrt(P2) = 0, as alpha is 0.")

    print(f"\nFinal optimal b = {b_optimal}")
    return b_optimal

# Example values. You can change these to test different scenarios.
# Scenario 1: Vertex is outside the interval
print("--- Scenario 1 ---")
solve_optimal_b(P1=5.0, P2=5.0, alpha=0.5)
print("\n" + "="*20 + "\n")

# Scenario 2: Vertex is inside the interval
print("--- Scenario 2 ---")
solve_optimal_b(P1=2.0, P2=8.0, alpha=0.5)
