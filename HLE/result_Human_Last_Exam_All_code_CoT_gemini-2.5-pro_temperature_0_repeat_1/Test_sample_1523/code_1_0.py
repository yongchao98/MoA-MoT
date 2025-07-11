import sympy

def get_inner_product_expression():
    """
    This function symbolically constructs the inner product (φ, Dφ)
    for a neutral scalar field in finite-temperature field theory.
    """
    # 1. Define the necessary mathematical symbols.
    # Coordinates: τ (imaginary time), x, y, z (space)
    # Parameters: m (mass), β (inverse temperature)
    tau, x, y, z = sympy.symbols('τ x y z', real=True)
    m, beta = sympy.symbols('m β', real=True, positive=True)

    # The field φ is a function of the coordinates.
    phi = sympy.Function('φ')(tau, x, y, z)

    # 2. Define the differential operator D = (-∂τ² - ∇² + m²)
    # Second derivative with respect to imaginary time τ
    d2_phi_dtau2 = phi.diff(tau, 2)

    # Laplacian operator ∇² applied to φ
    laplacian_phi = phi.diff(x, 2) + phi.diff(y, 2) + phi.diff(z, 2)

    # The operator D applied to the field φ
    D_phi = -d2_phi_dtau2 - laplacian_phi + m**2 * phi

    # 3. Construct the integrand for the inner product.
    # The integrand is φ multiplied by Dφ.
    integrand = phi * D_phi

    # 4. Represent the full inner product as a formal integral.
    # The integral is over imaginary time from 0 to β and all of 3D space.
    # We use sympy.Integral to represent this.
    # The spatial integral is over R^3, from -∞ to +∞ for x, y, z.
    inner_product = sympy.Integral(integrand, 
                                   (tau, 0, beta), 
                                   (x, -sympy.oo, sympy.oo), 
                                   (y, -sympy.oo, sympy.oo), 
                                   (z, -sympy.oo, sympy.oo))

    # 5. Print the final expression.
    # We use sympy's pretty printer for a clean mathematical output.
    print("In the functional integral representation for a neutral scalar field,")
    print("the operator D is the Euclidean Klein-Gordon operator: D = (-∂τ² - ∇² + m²).")
    print("\nThe inner product (φ, Dφ) is given by the integral:")
    
    # The final equation is built from the components defined above.
    # The expression for D_phi shows each term explicitly.
    # For example, -1 * ∂τ²(φ), -1 * ∇²(φ), and m**2 * φ.
    # The printout will render these components.
    sympy.pprint(inner_product, use_unicode=True)

if __name__ == '__main__':
    get_inner_product_expression()