def print_electromagnetic_solution():
    """
    This function prints the derived equations for the electric potential
    and electric field in the region outside the conducting sphere (r > R).
    The expressions are given in LaTeX format.
    """

    # Expression for the electric potential outside the sphere
    phi_out_str = (
        r"$\Phi(r, \theta) = "
        r"-E_0 \left( r - \frac{(\sigma_1 - \sigma_2) R^3}{(\sigma_1 + 2\sigma_2)r^2} \right) \cos\theta$"
    )

    # Expression for the electric field outside the sphere
    # This combines the radial (E_r) and polar (E_theta) components.
    e_out_str = (
        r"$\vec{E}(r, \theta) = "
        r"E_0 \left[ 1 + \frac{2(\sigma_1 - \sigma_2) R^3}{(\sigma_1 + 2\sigma_2) r^3} \right] \cos\theta \hat{r} "
        r"- E_0 \left[ 1 - \frac{(\sigma_1 - \sigma_2) R^3}{(\sigma_1 + 2\sigma_2) r^3} \right] \sin\theta \hat{\theta}$"
    )

    print("The electric potential for r > R is:")
    print(phi_out_str)
    print("\n" + "#" * 50 + "\n")
    print("The electric field for r > R is:")
    print(e_out_str)

# Execute the function to display the results.
print_electromagnetic_solution()