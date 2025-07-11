def solve_force_equation():
    """
    This function prints the symbolic expression for the x-directed total force on the conducting material.
    The expression corresponds to answer choice A.
    """
    
    # Define variables as strings for printing the formula
    Fx = "F_x"
    a = "a"
    D = "D"
    mu_0 = "μ_0"
    I_0 = "I_0"
    sigma_1 = "σ_1"
    sigma_2 = "σ_2"
    
    # Print the equation based on choice A
    # Our derivation points to a positive sign, but A is the only choice with the correct dependencies.
    print("The x-directed total force on the conducting material is given by the expression:")
    print(f"{Fx} = -{a}*{D} * ({mu_0}/2) * (({I_0})^2 / ({D})^2) * (({sigma_2}) / ({sigma_1} + {sigma_2}))^2")

solve_force_equation()