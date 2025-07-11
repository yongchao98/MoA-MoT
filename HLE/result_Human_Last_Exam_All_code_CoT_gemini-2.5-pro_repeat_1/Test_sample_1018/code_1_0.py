def compute_simplicial_volume_product():
    """
    Computes the simplicial volume of the product of two surfaces Sigma_g and Sigma_h.
    """
    # Genus of the two surfaces
    g = 31
    h = 17

    # Step 1 & 2: Calculate the simplicial volume for each surface.
    # The formula for a surface of genus g >= 2 is ||Sigma_g|| = 4g - 4.
    sv_g = 4 * g - 4
    sv_h = 4 * h - 4

    print(f"The simplicial volume of an oriented closed surface Sigma_g with genus g >= 2 is ||Sigma_g|| = 4*g - 4.")
    print("-" * 30)
    print(f"For the first surface, g = {g}:")
    print(f"||Sigma_{g}|| = 4 * {g} - 4 = {4*g} - 4 = {sv_g}")
    print("")
    print(f"For the second surface, h = {h}:")
    print(f"||Sigma_{h}|| = 4 * {h} - 4 = {4*h} - 4 = {sv_h}")
    print("-" * 30)

    # Step 3 & 4: Use the product formula for surfaces with genera g, h >= 2.
    # The formula is ||Sigma_g x Sigma_h|| = 6 * ||Sigma_g|| * ||Sigma_h||.
    sv_product = 6 * sv_g * sv_h

    print("The simplicial volume of the product of two surfaces Sigma_g x Sigma_h (for g,h >= 2) is:")
    print("||Sigma_g x Sigma_h|| = 6 * ||Sigma_g|| * ||Sigma_h||")
    print("")
    print("Substituting the calculated values into the formula:")
    print(f"||Sigma_{g} x Sigma_{h}|| = 6 * {sv_g} * {sv_h} = {sv_product}")

compute_simplicial_volume_product()