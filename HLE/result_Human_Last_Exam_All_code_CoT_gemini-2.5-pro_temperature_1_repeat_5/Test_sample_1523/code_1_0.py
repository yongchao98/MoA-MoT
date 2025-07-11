def print_inner_product_formula():
    """
    This function prints the mathematical formula for the inner product (ϕ, D_ϕ)
    in finite-temperature field theory for a neutral scalar field.
    """
    # The formula is constructed using unicode characters for mathematical symbols.
    # The inner product is an integral over the Euclidean spacetime volume.
    # At finite temperature T, the time integral is over the interval [0, β], where β = 1/T.
    # The spatial integral is over all of 3D space.
    # The integrand consists of the kinetic term (temporal and spatial derivatives squared)
    # and the mass term.
    # Note that the numbers 0, 3, and 2 are all present in the final equation.
    equation = "(ϕ, D_ϕ) = ∫₀^β dτ ∫ d³x [ (∂_τ ϕ)² + (∇ϕ)² + m²ϕ² ]"
    print(equation)

print_inner_product_formula()