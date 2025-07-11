def print_inner_product_expression():
    """
    This function prints the mathematical expression for the inner product (ϕ, D_ϕ)
    in finite-temperature field theory for a neutral scalar field.
    """

    # Define the components of the expression as strings
    integral_part = "∫[0 to β] dτ ∫ d³x"
    field_phi = "ϕ(τ, x)"
    operator_D_prefix = "["
    d2_dtau2_term = "-∂²/∂τ²"
    laplacian_term = "-∇²"
    mass_term = "m²"
    operator_D_suffix = "]"

    # The full expression for the inner product (ϕ, D_ϕ)
    # is the integral of ϕ multiplied by D applied to ϕ.
    print("In finite-temperature field theory, the inner product (ϕ, D_ϕ) is defined as:")
    print("")

    # We print each part of the equation to be explicit.
    # Final equation form: Integral * phi * (D_operator) * phi
    final_equation = (
        f"{integral_part} {field_phi} * "
        f"({operator_D_prefix} {d2_dtau2_term} {laplacian_term} + {mass_term} {operator_D_suffix}) "
        f"{field_phi}"
    )

    print(final_equation)
    print("")
    print("Where:")
    print("  - τ is the imaginary time, integrated from 0 to β = 1/T (T is temperature).")
    print("  - x represents the 3 spatial dimensions.")
    print("  - ϕ(τ, x) is the neutral scalar field.")
    print("  - ∂²/∂τ² is the second partial derivative with respect to imaginary time.")
    print("  - ∇² is the spatial Laplacian operator.")
    print("  - m is the mass of the scalar field.")

if __name__ == "__main__":
    print_inner_product_expression()