import math

def compute_simplicial_volume_product_surfaces():
    """
    Computes the simplicial volume of the product of two closed oriented surfaces.

    The formula for the simplicial volume of Sigma_g x Sigma_h for g, h >= 1 is:
    ||Sigma_g x Sigma_h|| = 6 * (g - 1) * (h - 1)
    """
    # Genus of the first surface
    g = 31
    # Genus of the second surface
    h = 17

    # Check if the genera are valid for the formula
    if g < 1 or h < 1:
        print("The formula is valid for genera g, h >= 1.")
        return

    # Calculate the components of the formula
    g_minus_1 = g - 1
    h_minus_1 = h - 1

    # Calculate the final simplicial volume
    simplicial_volume = 6 * g_minus_1 * h_minus_1

    # Print the full equation with the numerical values
    print(f"To compute the simplicial volume of Sigma_{g} x Sigma_{h}, we use the formula:")
    print(f"||Sigma_g x Sigma_h|| = 6 * (g - 1) * (h - 1)")
    print(f"Substituting g = {g} and h = {h}:")
    print(f"Result = 6 * ({g} - 1) * ({h} - 1) = 6 * {g_minus_1} * {h_minus_1} = {simplicial_volume}")

# Execute the function to print the result
compute_simplicial_volume_product_surfaces()