import math

def solve_contracts():
    """
    Calculates the separating equilibrium contracts for Chelsea and the handymen.
    
    The problem is solved by setting up Chelsea's profit maximization subject to
    the Incentive Compatibility (IC) and Individual Rationality (IR) constraints.
    
    The key steps are:
    1. Assume the IR constraints are binding to express prices in terms of values.
    2. Substitute these into the profit function and IC constraints.
    3. Maximize the profit function with respect to v_L, subject to the IC_H constraint.
       This yields the optimal v_L.
    4. Use the IC_L constraint, assumed to be binding, to solve for v_H.
    5. Calculate the corresponding prices p_L and p_H.
    """
    
    # Step 1 & 2: Set up the problem.
    # Chelsea's profit: E[pi] = (5/6)*(v_H - p_H) + (1/6)*(v_L - p_L)
    # Binding IR constraints:
    # p_H = v_H + 30
    # p_L = v_L^2 + 20 - (1/3)*v_L
    # Substituting these into the profit function gives:
    # E[pi] = (5/6)*(-30) + (1/6)*(v_L - (v_L^2 + 20 - v_L/3))
    # E[pi] = -25 + (1/6)*(-v_L^2 + (4/3)*v_L - 20)
    # To maximize E[pi], we must maximize f(v_L) = -v_L^2 + (4/3)*v_L.
    
    # Step 3: Solve for v_L.
    # The maximization is subject to the IC_H constraint: p_H - v_H >= p_L - v_L
    # -> 30 >= (v_L^2 + 20 - v_L/3) - v_L
    # -> 10 >= v_L^2 - (4/3)*v_L
    # -> v_L^2 - (4/3)*v_L - 10 <= 0
    # The unconstrained maximum of f(v_L) is found by setting its derivative to zero:
    # f'(v_L) = -2*v_L + 4/3 = 0  => v_L = 2/3.
    # We check if this solution satisfies the constraint:
    # (2/3)^2 - (4/3)*(2/3) - 10 = 4/9 - 8/9 - 10 = -4/9 - 10, which is <= 0.
    # So, the constraint is satisfied, and the optimal v_L is 2/3.
    v_L = 2 / 3
    
    # Calculate p_L from the binding IR_L constraint
    p_L = v_L**2 + 20 - (1/3) * v_L
    
    # Step 4: Solve for v_H.
    # We use the IC_L constraint: p_L - v_L^2 >= p_H - v_H^2
    # -> (v_L^2 + 20 - v_L/3) - v_L^2 >= (v_H + 30) - v_H^2
    # -> 20 - v_L/3 >= -v_H^2 + v_H + 30
    # -> v_H^2 - v_H - 10 - v_L/3 >= 0
    # Since profit does not depend on v_H, we assume this constraint is binding.
    # v_H^2 - v_H - (10 + v_L/3) = 0
    # v_H^2 - v_H - (10 + (2/3)/3) = 0
    # v_H^2 - v_H - (10 + 2/9) = 0
    # v_H^2 - v_H - 92/9 = 0
    # Using the quadratic formula for v_H = (-b + sqrt(b^2 - 4ac)) / 2a
    # a=1, b=-1, c=-92/9
    v_H_discriminant = (-1)**2 - 4 * 1 * (-92/9)
    # We take the positive root since value added should be positive.
    v_H = (1 + math.sqrt(v_H_discriminant)) / 2
    
    # Step 5: Calculate p_H from the binding IR_H constraint.
    p_H = v_H + 30
    
    # Output the results
    print("The separating equilibrium is defined by the following pair of contracts:")
    print("\nContract for Low Type (L):")
    print(f"(v_L, p_L) = ({v_L}, {p_L})")
    print(f"Approximately: (v_L, p_L) = ({v_L:.4f}, {p_L:.4f})")
    
    print("\nContract for High Type (H):")
    print(f"(v_H, p_H) = ({v_H}, {p_H})")
    print(f"Approximately: (v_H, p_H) = ({v_H:.4f}, {p_H:.4f})")

solve_contracts()