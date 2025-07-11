def solve_neron_severi_rank():
    """
    Calculates the smallest and largest possible rank of the Neron-Severi group
    of the 15th symmetric power of a genus 3 Riemann surface.
    """

    # The genus of the curve C
    g = 3

    # The degree of the symmetric power X = C^(d)
    d = 15

    # The rank of the Neron-Severi group of X is given by the formula:
    # rho(X) = rho(J(C)) + 1
    # where J(C) is the Jacobian of C.
    # We need to find the min and max possible values for rho(J(C)).

    # Smallest Rank Calculation:
    # For a generic curve C, the rank of the Neron-Severi group of its Jacobian
    # is minimal.
    # min rho(J(C)) = 1
    rho_J_min = 1
    rho_X_min = rho_J_min + 1
    print("Finding the smallest possible rank:")
    print(f"The minimum rank for the Jacobian's Neron-Severi group is rho(J(C)) = {rho_J_min}.")
    print(f"The smallest rank for X = C^(15) is rho(J(C)) + 1 = {rho_J_min} + 1 = {rho_X_min}.")
    print("-" * 20)

    # Largest Rank Calculation:
    # The rank of the Neron-Severi group of a g-dimensional abelian variety is
    # bounded by the Hodge number h^{1,1} = g*g. This maximum is achieved
    # for Jacobians with Complex Multiplication.
    # max rho(J(C)) = g^2
    rho_J_max = g * g
    rho_X_max = rho_J_max + 1
    print("Finding the largest possible rank:")
    print(f"The maximum rank for the Jacobian's Neron-Severi group is g*g = {g} * {g} = {rho_J_max}.")
    print(f"The largest rank for X = C^(15) is rho(J(C)) + 1 = {rho_J_max} + 1 = {rho_X_max}.")

solve_neron_severi_rank()