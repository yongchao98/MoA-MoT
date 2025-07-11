import sympy

def get_scalar_field_inner_product():
    """
    This function uses sympy to symbolically derive the integrand of the
    inner product (phi, D_phi) for a neutral scalar field in
    finite-temperature field theory.
    """
    # 1. Define the coordinates for 4D Euclidean spacetime.
    # tau is the imaginary time coordinate, x, y, z are spatial coordinates.
    tau, x, y, z = sympy.symbols('tau x y z', real=True)

    # 2. Define the neutral scalar field 'phi' as a function of these coordinates
    # and its mass 'm'.
    phi = sympy.Function('phi')(tau, x, y, z)
    m = sympy.Symbol('m', positive=True)

    # 3. The Euclidean action S_E can be written as S_E = 1/2 * (phi, D_phi),
    # where D is the Euclidean Klein-Gordon operator: D = (-box + m^2).
    # 'box' is the 4D Laplacian operator (d'Alembertian in Euclidean space).
    # box*phi = (d^2/dtau^2 + d^2/dx^2 + d^2/dy^2 + d^2/dz^2) * phi

    # Let's construct the term 'box*phi'.
    d2phi_dtau2 = phi.diff(tau, 2)
    d2phi_dx2 = phi.diff(x, 2)
    d2phi_dy2 = phi.diff(y, 2)
    d2phi_dz2 = phi.diff(z, 2)
    box_phi = d2phi_dtau2 + d2phi_dx2 + d2phi_dy2 + d2phi_dz2

    # 4. Now, construct the term D*phi.
    # The operator D acts on phi as D_phi = (-box + m**2) * phi
    D_phi = -box_phi + m**2 * phi

    # 5. The inner product (phi, D_phi) is the integral of phi * (D_phi)
    # over the 4D Euclidean volume. We will print the integrand.
    # This shows the full structure of the inner product before integration.
    integrand = phi * D_phi
    
    # We explicitly output each term as requested.
    # The term with the mass squared is m**2 * phi**2.
    mass_term = m**2 * phi**2
    # The term with the derivatives is -phi * box_phi.
    kinetic_term = -phi * box_phi

    print("In finite-temperature field theory for a neutral scalar field phi with mass m,")
    print("the Euclidean action S_E can be written in the form S_E = 1/2 * (phi, D_phi),")
    print("where D is the operator (-box + m^2) and box is the 4D Laplacian.\n")
    print("The inner product (phi, D_phi) is the integral of the following expression over 4D Euclidean spacetime:")
    print("\nIntegrand = phi * D_phi = phi * (-box + m^2)*phi\n")

    # The final equation is Integrand = mass_term + kinetic_term
    # We print it out symbolically. The only explicit number here is the exponent 2.
    print("Symbolic representation of the integrand:")
    sympy.pprint(integrand, use_unicode=False)

    print("\nThis can also be expressed (before integration by parts) in terms of first derivatives:")
    print("(del_mu phi)*(del^mu phi) + m^2*phi^2")

if __name__ == '__main__':
    get_scalar_field_inner_product()
