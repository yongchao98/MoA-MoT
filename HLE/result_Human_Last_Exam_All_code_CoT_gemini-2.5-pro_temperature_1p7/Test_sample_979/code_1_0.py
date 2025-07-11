def print_magnetic_field_solution():
    """
    This function prints the final expression for the magnetic field H(r, theta),
    deconstructing it into its constituent parts as requested.
    """

    # Components for the field inside the sphere (0 < r < R)
    H_in_term1 = "\\frac{2 \\mu_0}{\\mu}"
    H_in_term2 = "\\frac{K_0}{1 + \\frac{2 \\mu_0}{\\mu}}"
    H_in_vector = "\\hat{z}"
    H_in_full_expression = f"{H_in_term1} {H_in_term2} {H_in_vector}"

    # Components for the field outside the sphere (R < r < infinity)
    H_out_coeff = "\\frac{K_0}{1 + \\frac{2 \\mu_0}{\\mu}}"
    H_out_spatial = "\\frac{R^3}{r^3}"
    H_out_vector = "(2 \\cos \\theta \\hat{r} + \\sin \\theta \\hat{\theta})"
    H_out_full_expression = f"{H_out_coeff} {H_out_spatial} {H_out_vector}"
    
    print("The final derived magnetic field H(r, theta) is expressed in LaTeX format below:")
    print(r"$$")
    print(r"\vec{H}(r, \theta) = \begin{cases}")
    print(f"    {H_in_full_expression} & 0 < r < R \\\\")
    print(f"    {H_out_full_expression} & R < r < \\infty")
    print(r"\end{cases}")
    print(r"$$")

    print("\n--- Components for the field inside the sphere (0 < r < R) ---")
    print(f"Term 1 of coefficient: {H_in_term1}")
    print(f"Term 2 of coefficient: {H_in_term2}")
    print(f"Vector component: {H_in_vector}\n")

    print("--- Components for the field outside the sphere (R < r < infinity) ---")
    print(f"Coefficient part: {H_out_coeff}")
    print(f"Spatial dependence: {H_out_spatial}")
    print(f"Angular vector part: {H_out_vector}")

# Execute the function to display the results
print_magnetic_field_solution()