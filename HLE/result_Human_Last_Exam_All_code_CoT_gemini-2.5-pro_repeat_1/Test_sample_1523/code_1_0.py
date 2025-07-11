import sympy

def solve_inner_product():
    """
    Symbolically computes and displays the inner product (ϕ, D_ϕ) for a neutral
    scalar field in finite-temperature field theory.
    """
    # Define spacetime coordinates and mass as symbolic variables
    # τ is imaginary time, x_1, x_2, x_3 are spatial coordinates
    tau, x1, x2, x3 = sympy.symbols('τ x_1 x_2 x_3')
    m = sympy.Symbol('m', positive=True)

    # Define the scalar field ϕ as a function of the spacetime coordinates
    phi = sympy.Function('ϕ')(tau, x1, x2, x3)

    # Define the 4-dimensional Euclidean Laplacian operator acting on ϕ.
    # The Laplacian ∂² is (∂²/∂τ² + ∂²/∂x₁² + ∂²/∂x₂² + ∂²/∂x₃²)
    laplacian_phi = (
        sympy.Derivative(phi, tau, 2) +
        sympy.Derivative(phi, x1, 2) +
        sympy.Derivative(phi, x2, 2) +
        sympy.Derivative(phi, x3, 2)
    )

    # The operator D in the quadratic action S₀ = ½ ∫ ϕ D ϕ is D = [-∂² + m²].
    # We now find the expression for D applied to ϕ.
    D_phi = -laplacian_phi + m**2 * phi

    # The inner product (ϕ, D_ϕ) is defined as the integral of ϕ * (D_ϕ)
    # over 4D Euclidean spacetime.
    integrand = phi * D_phi

    # We represent the full inner product symbolically using sympy.Integral.
    # The integral is over d⁴x which is represented by the coordinates.
    inner_product_expression = sympy.Integral(integrand, (tau, x1, x2, x3))

    # Print the final symbolic expression for the inner product.
    # sympy.pprint provides a nicely formatted output.
    print("The symbolic expression for the inner product (ϕ, D_ϕ) is:")
    sympy.pprint(inner_product_expression)

    # To fulfill the requirement to "output each number in the final equation",
    # we can explicitly print the terms of the integrand. The numbers here are the
    # exponents '2' in the mass term and the order '2' of the derivatives.
    print("\nThe integrand contains the following terms:")
    print(f"Mass term component: ϕ * m**2 * ϕ")
    print(f"Kinetic term component (from time derivative): -ϕ * d**2/dτ**2(ϕ)")
    print(f"Kinetic term component (from spatial derivatives): -ϕ * ∇**2(ϕ)")

solve_inner_product()