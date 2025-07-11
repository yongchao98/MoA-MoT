def solve_geometry_problem():
    """
    This function explains the derivation and prints the final equation.
    """
    
    # Define string representations for the geometric segments
    MA = "MA"
    MB = "MB"
    MG = "MG"
    MH = "MH"
    
    # Based on the geometric proof, we found that G and H are the same point.
    # Therefore, the lengths MG and MH are equal.
    # The result of MG - MH is 0.
    result = 0
    
    # The problem asks to express the result in terms of MA and MB.
    # We can write 0 as 0 * MA + 0 * MB, which evaluates to 0.
    # The final equation will show all the numbers involved.
    
    coeff_ma = 0
    coeff_mb = 0
    
    print("The derivation shows that G and H are the same point, so MG = MH.")
    print("Therefore, the expression MG - MH equals 0.")
    print("The final equation, showing all numbers as requested, is:")
    
    # Print the equation with coefficients for MA and MB
    # This satisfies the requirement to "output each number in the final equation".
    print(f"{MG} - {MH} = {coeff_ma} * {MA} + {coeff_mb} * {MB}")

    # Print the simplified final result
    print(f"Which simplifies to: {MG} - {MH} = {result}")

solve_geometry_problem()
