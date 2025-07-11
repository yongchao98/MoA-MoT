def display_scalar_inner_product():
    """
    This function explains and displays the inner product (ϕ, D_ϕ)
    for a neutral scalar field in finite-temperature field theory.
    """
    print("In finite-temperature field theory, the kinetic part of the Euclidean action S_E for a scalar field ϕ can be written in the form S_E = ½(ϕ, Dϕ).")
    print("This allows us to identify the differential operator D and the associated inner product.\n")
    print("The operator D is found to be: D = -∂_τ² - ∇² + m²\n")
    print("The inner product (f, g) is defined as ∫ dτ d³x f(x)g(x).")
    print("Therefore, the inner product (ϕ, Dϕ) is given by the following expression.")
    print("The equation explicitly shows each numerical coefficient as requested.\n")

    # Define symbolic components for clarity in the final output
    phi_symbol = "ϕ"
    dtau2_phi_symbol = "∂_τ²ϕ"
    nabla2_phi_symbol = "∇²ϕ"
    m_symbol = "m"

    # Coefficients of the operator D
    coeff_dtau2 = -1
    coeff_nabla2 = -1
    coeff_m2 = 1

    # Construct the string for the integrand, emphasizing the coefficients.
    integrand = (
        f"{phi_symbol} * ("
        f"({coeff_dtau2}) * {dtau2_phi_symbol} + "
        f"({coeff_nabla2}) * {nabla2_phi_symbol} + "
        f"({coeff_m2}) * {m_symbol}² * {phi_symbol}"
        ")"
    )

    # Construct the final string for the full inner product.
    final_expression = f"(ϕ, Dϕ) = ∫ dτ d³x [ {integrand} ]"

    print("="*60)
    print(final_expression)
    print("="*60)

if __name__ == '__main__':
    display_scalar_inner_product()