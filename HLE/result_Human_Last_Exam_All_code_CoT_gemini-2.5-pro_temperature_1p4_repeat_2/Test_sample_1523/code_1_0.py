import sympy

def solve_task():
    """
    This function uses the sympy library to symbolically represent the inner product
    (ϕ, Dϕ) in finite-temperature field theory for a neutral scalar field.
    """

    # Set up pretty printing for a clear mathematical output
    sympy.init_printing(use_unicode=True)

    # --- 1. Define the mathematical objects ---

    # Define Euclidean spacetime coordinates: τ (imaginary time) and x, y, z (space).
    tau, x, y, z = sympy.symbols('τ x y z', real=True)

    # The neutral scalar field ϕ is a function of spacetime coordinates.
    phi = sympy.Function('ϕ')(x, y, z, tau)

    # Define the parameters of the theory: m (mass) and β (inverse temperature).
    m = sympy.Symbol('m', real=True, positive=True)
    beta = sympy.Symbol('β', real=True, positive=True)


    # --- 2. Define the operator D acting on the field ϕ ---

    # The operator D is the Euclidean Klein-Gordon operator: D = -∂²/∂τ² - ∇² + m².
    # We construct its action on the field ϕ.

    # Second derivative with respect to imaginary time τ
    d2_dtau2_phi = sympy.diff(phi, tau, 2)

    # Laplacian in spatial coordinates (x, y, z)
    nabla_squared_phi = sympy.laplacian(phi, (x, y, z))

    # The expression for D applied to ϕ (D_ϕ)
    D_phi = -d2_dtau2_phi - nabla_squared_phi + sympy.S(1) * m**2 * phi


    # --- 3. Define the inner product (ϕ, D_ϕ) ---

    # The inner product is the integral of ϕ * (D_ϕ) over Euclidean spacetime.
    # The integration domain is over all of 3D space and from 0 to β for imaginary time.

    # The integrand of the inner product
    integrand = phi * D_phi

    # The full inner product expression as a Sympy Integral object.
    # The numbers in the equation, such as the coefficient -1 for the derivatives
    # and +1 for the mass term, are explicitly represented.
    inner_product = sympy.Integral(integrand,
                                  (tau, 0, beta),
                                  (x, -sympy.oo, sympy.oo),
                                  (y, -sympy.oo, sympy.oo),
                                  (z, -sympy.oo, sympy.oo))


    # --- 4. Print the final expression ---

    print("The inner product (ϕ, Dϕ) for a neutral scalar field is given by the expression:")
    print("\n(ϕ, Dϕ) = ")
    sympy.pprint(inner_product)

solve_task()