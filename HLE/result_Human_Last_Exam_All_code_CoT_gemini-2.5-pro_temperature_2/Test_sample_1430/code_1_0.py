def solve_neron_severi_rank():
    """
    Calculates the smallest and largest possible rank of the Neron-Severi group
    for the 15th symmetric power of a genus 3 Riemann surface.
    """
    
    # The genus of the Riemann surface C
    g = 3
    
    print(f"The curve C has genus g = {g}.")
    print("The rank of the Neron-Severi group of X = C^(15) is given by the formula: rho(X) = rho(Jac(C)) + 1.")
    print("The rank of the Neron-Severi group of the Jacobian, rho(Jac(C)), is bounded by 1 <= rho(Jac(C)) <= g^2.")
    print("-" * 30)

    # Smallest possible rank
    # This occurs for a generic curve C, where rho(Jac(C)) = 1.
    rho_jac_min = 1
    rho_X_min = rho_jac_min + 1
    print("Calculating the smallest possible rank:")
    print(f"The minimum rank for rho(Jac(C)) is {rho_jac_min}.")
    print(f"Smallest rank of NS(X) = {rho_jac_min} + 1 = {rho_X_min}")
    print("-" * 30)
    
    # Largest possible rank
    # This occurs for a special curve C, where rho(Jac(C)) = g^2.
    rho_jac_max = g**2
    rho_X_max = rho_jac_max + 1
    print("Calculating the largest possible rank:")
    print(f"The maximum rank for rho(Jac(C)) is g^2 = {g}^2 = {rho_jac_max}.")
    print(f"Largest rank of NS(X) = {rho_jac_max} + 1 = {rho_X_max}")
    print("-" * 30)

solve_neron_severi_rank()
