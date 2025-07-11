import math

def find_optimal_b(alpha, P1, P2):
    """
    Calculates the optimal feedback adjustment factor b.

    Args:
        alpha (float): The weather-induced correlation between noise measurements.
        P1 (float): Power (as a dimensionless SNR) used in the first transmission.
        P2 (float): Power (as a dimensionless SNR) used in the second transmission.

    Returns:
        float: The optimal feedback adjustment factor b.
    """
    # The unconstrained optimal value for b is derived as -alpha * P1
    b_unconstrained = -alpha * P1

    # The value of b is constrained by the power of the second transmission: b^2 <= P2
    if b_unconstrained**2 <= P2:
        # The unconstrained optimum is feasible
        b_optimal = b_unconstrained
        print("The optimal feedback factor b is within the power constraint.")
        print(f"b = -alpha * P1 = -({alpha}) * ({P1}) = {b_optimal}")
    else:
        # The unconstrained optimum is not feasible; the optimum lies on the boundary.
        # The boundary point is determined by the sign of -alpha * P1, which is -sgn(alpha)
        # sgn(alpha) is math.copysign(1, alpha) if alpha is not 0
        if alpha == 0:
             # If alpha is 0, b_unconstrained is 0, which is always in [-sqrt(P2), sqrt(P2)]
             # This case is handled above, but included for completeness.
            sgn_alpha = 0
        else:
            sgn_alpha = math.copysign(1, alpha)
        
        b_optimal = -sgn_alpha * math.sqrt(P2)
        print("The unconstrained optimum is outside the power constraint.")
        print(f"The optimal b is on the boundary: b = -sgn(alpha) * sqrt(P2) = -sgn({alpha}) * sqrt({P2}) = {b_optimal}")
        
    return b_optimal

if __name__ == '__main__':
    # Example 1: Unconstrained case
    print("--- Example 1: Unconstrained Optimum is Feasible ---")
    alpha1 = 0.5
    P1_1 = 1.0
    P2_1 = 1.5
    print(f"Given parameters: alpha = {alpha1}, P1 = {P1_1}, P2 = {P2_1}")
    print(f"Constraint check: (-alpha * P1)^2 <= P2  =>  ({-alpha1*P1_1})^2 <= {P2_1}  =>  {(-alpha1*P1_1)**2} <= {P2_1} -> {(-alpha1*P1_1)**2 <= P2_1}")
    b_opt1 = find_optimal_b(alpha1, P1_1, P2_1)
    print(f"\nFinal optimal b = {b_opt1}")
    
    print("\n" + "="*50 + "\n")

    # Example 2: Constrained (boundary) case
    print("--- Example 2: Unconstrained Optimum is Not Feasible ---")
    alpha2 = 0.8
    P1_2 = 2.0
    P2_2 = 1.0
    print(f"Given parameters: alpha = {alpha2}, P1 = {P1_2}, P2 = {P2_2}")
    print(f"Constraint check: (-alpha * P1)^2 <= P2  =>  ({-alpha2*P1_2})^2 <= {P2_2}  =>  {round((-alpha2*P1_2)**2, 2)} <= {P2_2} -> {(-alpha2*P1_2)**2 <= P2_2}")
    b_opt2 = find_optimal_b(alpha2, P1_2, P2_2)
    print(f"\nFinal optimal b = {b_opt2}")
