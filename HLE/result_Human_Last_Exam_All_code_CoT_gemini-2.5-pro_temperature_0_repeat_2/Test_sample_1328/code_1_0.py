import math

def find_optimal_b(P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    Args:
        P1 (float): Power used in the first transmission.
        P2 (float): Power used in the second transmission.
        alpha (float): Weather-induced correlation coefficient.
    """
    print(f"Given parameters:")
    print(f"P1 = {P1}")
    print(f"P2 = {P2}")
    print(f"alpha = {alpha}\n")

    # The optimization problem is to maximize f(b) = -b^2 - 2*alpha*P1*b + C
    # subject to the constraint b^2 <= P2.

    # The unconstrained maximum is at b = -alpha * P1.
    # We check if this solution is within the feasible region defined by the constraint.
    
    condition_value = alpha**2 * P1**2
    
    print(f"Checking condition: P2 >= alpha^2 * P1^2")
    print(f"Is {P2} >= {alpha**2} * {P1**2}? -> Is {P2} >= {condition_value}?")

    if P2 >= condition_value:
        # The unconstrained optimum is feasible.
        optimal_b = -alpha * P1
        print("\nCondition is TRUE.")
        print("The optimal feedback factor 'b' is given by the equation: b = -alpha * P1")
        print(f"b = -({alpha}) * {P1} = {optimal_b}")
    else:
        # The unconstrained optimum is outside the feasible region.
        # The optimum is at the boundary of the constraint b^2 <= P2.
        # The sign of b is chosen to be opposite to the sign of alpha to maximize the objective.
        optimal_b = -math.copysign(1.0, alpha) * math.sqrt(P2)
        print("\nCondition is FALSE.")
        print("The optimal feedback factor 'b' is limited by P2 and is given by: b = -sign(alpha) * sqrt(P2)")
        print(f"b = -sign({alpha}) * sqrt({P2}) = {optimal_b}")
        
    return optimal_b

# --- Example Cases ---

# Case 1: The unconstrained solution is feasible.
print("--- Case 1: Unconstrained Solution is Feasible ---")
find_optimal_b(P1=2.0, P2=5.0, alpha=0.5)

print("\n" + "="*50 + "\n")

# Case 2: The unconstrained solution is not feasible, optimum is on the boundary.
print("--- Case 2: Unconstrained Solution is Not Feasible ---")
find_optimal_b(P1=10.0, P2=5.0, alpha=0.5)
