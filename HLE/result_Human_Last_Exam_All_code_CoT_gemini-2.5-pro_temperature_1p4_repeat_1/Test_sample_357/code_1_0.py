import math

def solve_handyman_contracts():
    """
    This script calculates the separating equilibrium contracts for Chelsea's handyman problem.
    """
    # Step 1: Determine v_L.
    # Chelsea's profit maximization leads to choosing v_L to maximize
    # -v_L^2 + (4/3)*v_L - 20.
    # The maximum is at v_L = (4/3) / 2 = 2/3.
    v_L = 2/3

    # Step 2: Determine p_L from the binding Individual Rationality (IR) constraint for the Low type.
    # p_L - v_L^2 = 20 - (1/3)*v_L
    p_L = v_L**2 + 20 - (1/3) * v_L

    # Step 3: Determine v_H.
    # This is determined from the binding Incentive Compatibility (IC) constraint for the Low type
    # and the binding IR constraint for the High type.
    # The resulting quadratic equation for v_H is: v_H^2 - v_H - (10 + (1/3)*v_L) = 0.
    # We solve for v_H using the quadratic formula, taking the positive root.
    a = 1
    b = -1
    c = -(10 + (1/3) * v_L)
    discriminant = b**2 - 4 * a * c
    v_H = (-b + math.sqrt(discriminant)) / (2 * a)

    # Step 4: Determine p_H from the binding IR constraint for the High type.
    # p_H - v_H = 30
    p_H = v_H + 30

    # Output the results
    print("The separating equilibrium consists of two contracts:")
    print("\n1. Contract for the Low-type handyman (L):")
    print(f"   Value added (v_L): {v_L:.4f}")
    print(f"   Price (p_L): {p_L:.4f}")

    print("\n2. Contract for the High-type handyman (H):")
    print(f"   Value added (v_H): {v_H:.4f}")
    print(f"   Price (p_H): {p_H:.4f}")
    
    # Return the values for the final answer block.
    # For symbolic answer:
    # v_L = 2/3
    # p_L = 182/9
    # v_H = (3 + sqrt(377))/6
    # p_H = (183 + sqrt(377))/6
    return (v_L, p_L, v_H, p_H)

# Execute the function
solve_handyman_contracts()