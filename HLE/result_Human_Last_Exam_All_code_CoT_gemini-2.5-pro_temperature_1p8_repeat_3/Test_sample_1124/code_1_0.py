def modified_logistic_map():
    """
    This function demonstrates a modified logistic map with a specific equilibrium point.
    """
    # System parameters
    R = 3.57
    C = 4/3

    # The theoretical equilibrium point X_eq = C - 1/R
    x_eq = C - 1/R

    print(f"The modified logistic map is: X_n+1 = R * X_n * (C - X_n)")
    print(f"With R = {R} and C = 4/3, the theoretical equilibrium point is {x_eq:.5f}")
    print("\nDemonstrating the equilibrium property: X_eq = f(X_eq, R)")
    
    # Start with the calculated equilibrium point
    x_n = x_eq
    
    # Calculate one iteration to show it's a fixed point
    x_n_plus_1 = R * x_n * (C - x_n)
    
    # Display the numbers in the final equation at equilibrium
    print("--------------------------------------------------")
    print("Final equation with calculated equilibrium value:")
    print(f"Input X_n: {x_n:.5f}")
    
    # Breaking down the calculation for clarity
    term1 = R
    term2 = x_n
    term3_calc = C - x_n
    term_mult = term1 * term2
    
    print(f"{x_n_plus_1:.5f} = {term1} * {term2:.5f} * ({C:.5f} - {term2:.5f})")
    print(f"{x_n_plus_1:.5f} = {term1} * {term2:.5f} * {term3_calc:.5f}")
    print(f"{x_n_plus_1:.5f} = {term_mult:.5f} * {term3_calc:.5f}")
    print(f"Output X_n+1: {x_n_plus_1:.5f}")
    print("--------------------------------------------------")

    # Final check
    if abs(x_n - x_n_plus_1) < 1e-9:
        print("As you can see, X_n+1 is approximately equal to X_n.")
        print("This confirms that the point is an equilibrium point of the modified map.")
    else:
        print("There was an error in the calculation, the point is not at equilibrium.")


modified_logistic_map()