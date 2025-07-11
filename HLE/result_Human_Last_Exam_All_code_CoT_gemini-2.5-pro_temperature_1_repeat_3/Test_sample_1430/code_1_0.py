import math

def solve_neron_severi_rank():
    """
    Calculates the smallest and largest possible rank of the Neron-Severi group
    of the 15th symmetric power of a genus 3 Riemann surface.
    """
    # The genus of the Riemann surface C
    g = 3
    
    # The degree of the symmetric power X = C^(d)
    d = 15

    print(f"The problem is to find the min/max rank of the Neron-Severi group of X = C^({d}) for a curve C of genus g = {g}.")
    print("The rank is given by the formula: rho(X) = rho(J(C)) + 1, where J(C) is the Jacobian of C.")
    print("-" * 20)
    
    # Smallest possible rank of NS(J(C)) for a generic curve of genus g
    # This corresponds to End(J(C)) = Z.
    min_rho_J = 1
    
    # Largest possible rank of NS(J(C)) for a special curve of genus g
    # This corresponds to a Jacobian with maximal complex multiplication. The rank is bounded by g^2.
    max_rho_J = g**2

    print(f"For a genus g = {g} curve, the rank of the Neron-Severi group of its Jacobian, rho(J(C)), is bounded:")
    print(f"Smallest rho(J(C)) = {min_rho_J} (for a generic curve).")
    print(f"Largest rho(J(C)) = g^2 = {g}^2 = {max_rho_J} (for a special CM curve).")
    print("-" * 20)

    # Calculate the smallest rank for X
    min_rho_X = min_rho_J + 1

    # Calculate the largest rank for X
    max_rho_X = max_rho_J + 1
    
    print("Now we can calculate the smallest and largest ranks for X = C^(15):")
    print(f"Smallest rank of NS(X) = {min_rho_J} + 1 = {min_rho_X}")
    print(f"Largest rank of NS(X) = {max_rho_J} + 1 = {max_rho_X}")

solve_neron_severi_rank()