def compute_simplicial_volume_product():
    """
    Computes the simplicial volume of the product of two closed oriented surfaces.
    """
    # Define the genus for each surface
    g = 31
    h = 17

    print("Step 1: Compute the simplicial volume of the first surface, Sigma_g.")
    print("The formula is ||Sigma_g|| = 4g - 4.")
    # Calculate the simplicial volume of Sigma_g
    sv_g = 4 * g - 4
    print(f"For g = {g}, ||Sigma_{g}|| = 4 * {g} - 4 = {4*g} - 4 = {sv_g}")
    print("-" * 30)

    print("Step 2: Compute the simplicial volume of the second surface, Sigma_h.")
    print("The formula is ||Sigma_h|| = 4h - 4.")
    # Calculate the simplicial volume of Sigma_h
    sv_h = 4 * h - 4
    print(f"For h = {h}, ||Sigma_{h}|| = 4 * {h} - 4 = {4*h} - 4 = {sv_h}")
    print("-" * 30)

    print("Step 3: Compute the simplicial volume of the product, Sigma_g x Sigma_h.")
    print("The formula is ||Sigma_g x Sigma_h|| = 6 * ||Sigma_g|| * ||Sigma_h||.")
    # Calculate the simplicial volume of the product
    sv_product = 6 * sv_g * sv_h
    print(f"The final simplicial volume is:")
    print(f"||\Sigma_{{{g}}} x \Sigma_{{{h}}}|| = 6 * {sv_g} * {sv_h} = {sv_product}")


compute_simplicial_volume_product()