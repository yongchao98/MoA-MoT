import sympy as sp

def calculate_inner_product_integrand():
    """
    Symbolically calculates and prints the integrand for the inner product (ϕ, Dϕ)
    in scalar field theory.
    """
    # Define the coordinates and constants as symbolic variables
    tau, x, y, z = sp.symbols('tau x y z')
    m = sp.Symbol('m')

    # Define the scalar field phi as a function of the spacetime coordinates
    phi = sp.Function('phi')(tau, x, y, z)

    # Calculate the partial derivatives of the field with respect to each coordinate
    d_phi_d_tau = sp.diff(phi, tau)
    d_phi_d_x = sp.diff(phi, x)
    d_phi_d_y = sp.diff(phi, y)
    d_phi_d_z = sp.diff(phi, z)

    # Construct the squared kinetic term, (∂_μ ϕ)²
    # Note: (∂_μ ϕ)² = (∂_τ ϕ)² + (∂_x ϕ)² + (∂_y ϕ)² + (∂_z ϕ)²
    kinetic_term_sq = d_phi_d_tau**2 + d_phi_d_x**2 + d_phi_d_y**2 + d_phi_d_z**2

    # Construct the mass term, m² ϕ²
    mass_term = m**2 * phi**2

    # The full integrand is the sum of the kinetic and mass terms.
    integrand = kinetic_term_sq + mass_term

    # Print the result. The inner product is the integral of this expression
    # over Euclidean spacetime (dτ d³x).
    print("The inner product (ϕ, Dϕ) is the integral over Euclidean spacetime of the following expression:")
    print("Integrand = ")
    sp.pretty_print(integrand)
    print("\nwhere:")
    print("- D is the Euclidean Klein-Gordon operator, D = -∂_μ² + m².")
    print("- (∂(phi)/∂tau)**2 + (∂(phi)/∂x)**2 + ... represents (∂_μ ϕ)².")
    print("- The final inner product is equal to 2 * S_free, where S_free is the action for the free field.")


if __name__ == '__main__':
    calculate_inner_product_integrand()