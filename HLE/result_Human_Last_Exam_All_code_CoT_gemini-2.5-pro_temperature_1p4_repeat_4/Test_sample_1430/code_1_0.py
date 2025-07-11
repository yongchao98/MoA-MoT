import math

def solve_neron_severi_rank():
    """
    Calculates the smallest and largest possible rank of the Neron-Severi group
    of the 15th symmetric power of a genus 3 Riemann surface.
    """

    # --- Problem setup ---
    # C is a genus g Riemann surface
    g = 3
    # X is the d-th symmetric power of C
    d = 15

    print(f"The problem is to find the min and max rank of NS(X) for X = Sym^{d}(C), where C is a curve of genus g={g}.")
    print("-" * 30)

    # --- Key relationship ---
    # The rank of the Neron-Severi group of X, rho(X), is related to the
    # rank of the Neron-Severi group of the Jacobian of C, rho(Jac(C)).
    # The formula rho(Sym^d(C)) = rho(Jac(C)) + 1 holds if d >= 2g - 1.
    condition = 2 * g - 1
    if d >= condition:
        print(f"The condition d >= 2g - 1 is met: {d} >= {condition}.")
        print("We can use the formula: rank(NS(X)) = rank(NS(Jac(C))) + 1.")
    else:
        print(f"The condition d >= 2g - 1 is NOT met: {d} < {condition}. The method is not applicable.")
        return

    print("-" * 30)

    # --- Smallest Possible Rank ---
    # For a generic curve C, the rank of the Neron-Severi group of its Jacobian
    # is minimal.
    rho_J_min = 1

    smallest_rank = rho_J_min + 1
    
    print("Calculation for the SMALLEST rank:")
    print(f"The minimum possible rank of NS(Jac(C)) for a genus {g} curve is {rho_J_min}.")
    print(f"Smallest rank of NS(X) = rank(NS(Jac(C)))_min + 1 = {rho_J_min} + 1 = {smallest_rank}")
    
    print("-" * 30)

    # --- Largest Possible Rank ---
    # The rank of the Neron-Severi group of an abelian variety of dimension g
    # is bounded by g^2. This maximum is achieved for special Jacobians.
    rho_J_max = g**2
    
    largest_rank = rho_J_max + 1

    print("Calculation for the LARGEST rank:")
    print(f"The maximum possible rank of NS(Jac(C)) for a genus {g} curve is g^2 = {g}^2 = {rho_J_max}.")
    print(f"Largest rank of NS(X) = rank(NS(Jac(C)))_max + 1 = {rho_J_max} + 1 = {largest_rank}")
    print("-" * 30)


# Execute the function to print the solution
solve_neron_severi_rank()