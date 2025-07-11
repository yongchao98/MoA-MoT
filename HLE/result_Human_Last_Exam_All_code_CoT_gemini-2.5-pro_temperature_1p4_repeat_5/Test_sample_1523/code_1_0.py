def print_inner_product_expression():
    """
    This function prints the final expression for the inner product
    (phi, D_phi) for a neutral scalar field in finite-temperature
    field theory.
    """

    expression = "(ϕ, Dϕ) = ∫d⁴x [ 1*(∂_τ ϕ)² + 1*(∇ϕ)² + 1*m²ϕ² ]"

    print("In the functional integral representation of the partition function, the action S[ϕ] is often written as S[ϕ] = 1/2 * (ϕ, Dϕ).")
    print("The inner product you asked for, (ϕ, Dϕ), is therefore equal to 2*S[ϕ].")
    print("\nIts explicit form is:")
    print(expression)
    print("\nWhere:")
    print("  - (ϕ, Dϕ): The inner product.")
    print("  - ∫d⁴x: The integral over Euclidean spacetime (from τ=0 to β, and over all 3D space).")
    print("  - (∂_τ ϕ)²: The square of the derivative of the field ϕ with respect to imaginary time τ.")
    print("  - (∇ϕ)²: The square of the spatial gradient of the field ϕ.")
    print("  - m: The mass of the scalar particle.")
    print("  - ϕ²: The square of the field ϕ.")
    print("  - The numbers '1' are the coefficients for each term in the integrand.")

if __name__ == "__main__":
    print_inner_product_expression()