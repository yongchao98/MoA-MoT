def solve_milp_constraints():
    """
    This function generates and prints the two additional MILP inequalities.
    """
    # Define the bounds for x. You can change these to any values.
    l = -10
    u = 10

    # The logic for choosing a and b is assumed to be handled by the MILP solver
    # based on the region of x:
    # a = 1 if x >= 1, a = 0 if x < 1
    # When a = 0: b = 1 if x >= 0, b = 0 if x < 0

    # The two new inequalities are derived using the big-M method to adjust the
    # lower bound on y when a = 0.

    # M1 is the coefficient for the first new inequality. Its value must be at
    # least max(0, -l) to ensure the constraint is redundant when not active.
    M1 = max(0, -l)

    # M2 is the coefficient for the second new inequality. Its value must be at
    # least 1.
    M2 = 1

    # The first inequality is derived from the condition that if a=0 and b=1, then y >= 0.
    # The general form is y >= 0 - M1*a - M1*(1-b).
    # After distributing, this is y >= -M1*a - M1 + M1*b.
    ineq1 = f"y >= {-M1}*a - {M1} + {M1}*b"

    # The second inequality is derived from the condition that if a=0 and b=0, then y >= x.
    # The general form is y >= x - M2*a - M2*b.
    ineq2 = f"y >= x - {M2}*a - {M2}*b"

    # Print the final inequalities in the requested format.
    print(f"{ineq1}, {ineq2}")

solve_milp_constraints()