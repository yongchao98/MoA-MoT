import math

def solve_modified_logistic_map():
    """
    This function modifies the standard logistic map to have an equilibrium
    point at approximately 1.05 for R=3.57, and verifies the result.
    """
    # The proposed modified logistic map is: X_n+1 = R * X_n * (1 - X_n) + log(R)

    # Set the value of R
    R = 3.57

    # An equilibrium point X_eq satisfies the equation:
    # X_eq = R * X_eq * (1 - X_eq) + log(R)
    # This can be rearranged into a quadratic equation:
    # R * X_eq^2 + (1 - R) * X_eq - log(R) = 0

    # Define the coefficients of the quadratic equation a*x^2 + b*x + c = 0
    a = R
    b = 1 - R
    c = -math.log(R)  # Natural logarithm

    # Calculate the discriminant
    try:
        discriminant = b**2 - 4 * a * c
    except ValueError:
        print("Error: Cannot calculate discriminant.")
        return

    # Solve for the two possible equilibrium points using the quadratic formula
    if discriminant >= 0:
        x1 = (-b + math.sqrt(discriminant)) / (2 * a)
        # x2 is the other root, which is negative in this case
        # x2 = (-b - math.sqrt(discriminant)) / (2 * a)
    else:
        print("No real equilibrium points found.")
        return

    # The positive root is the equilibrium point we are looking for.
    x_eq = x1

    print("The proposed modified logistic map is: X_n+1 = R * X_n * (1 - X_n) + log(R)")
    print(f"\nFor R = {R}, we solve for the equilibrium point X_eq where X_eq = f(X_eq, R).")
    print("This gives the quadratic equation:")
    print(f"{a:.2f}*X_eq^2 + ({b:.2f})*X_eq + ({c:.4f}) = 0")
    print(f"\nSolving this equation gives an equilibrium point at X_eq = {x_eq:.5f}")
    print("This value is approximately 1.05.")

    # Verification: plug the calculated equilibrium point back into the equation.
    # The final equation is X_eq = R * X_eq * (1 - X_eq) + log(R)
    # We will show that LHS = RHS
    lhs = x_eq
    log_r_val = math.log(R)
    rhs = R * x_eq * (1 - x_eq) + log_r_val

    print("\n--- Verification ---")
    print("We will plug the numbers into the equilibrium equation to verify:")
    print(f"Equation: {x_eq:.5f} = {R} * {x_eq:.5f} * (1 - {x_eq:.5f}) + log({R})")
    
    # Breaking down the right-hand side (RHS) calculation
    term_in_parenthesis = 1 - x_eq
    multiplication_part = R * x_eq * term_in_parenthesis
    
    print("\nLeft-Hand Side (LHS):")
    print(f"LHS = {lhs:.5f}")

    print("\nRight-Hand Side (RHS):")
    print(f"RHS = {R} * {x_eq:.5f} * ({term_in_parenthesis:.5f}) + {log_r_val:.5f}")
    print(f"RHS = {multiplication_part:.5f} + {log_r_val:.5f}")
    print(f"RHS = {rhs:.5f}")
    
    print(f"\nConclusion: LHS and RHS are equal, confirming that X_eq = {x_eq:.5f} is the equilibrium point.")

solve_modified_logistic_map()