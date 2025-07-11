import sympy as sp

def display_inner_product_integrand():
    """
    This function uses the sympy library to symbolically represent the integrand
    of the inner product (ϕ, D_ϕ) for a neutral scalar field.
    """
    
    # Define the symbols required for the expression
    # m is the mass of the field
    m = sp.Symbol('m', real=True)
    
    # τ represents imaginary time
    tau = sp.Symbol('τ', real=True)
    
    # ϕ is the scalar field, represented as a function. The '...' indicates
    # it is also a function of unspecified spatial coordinates.
    phi = sp.Function('ϕ')(tau, ...)
    
    # (∇ϕ)² is the square of the spatial gradient of the field.
    # We represent it as a single symbol for generality across any number of dimensions.
    grad_phi_sq = sp.Symbol('(∇ϕ)²')

    # Construct the three terms of the integrand based on the derivation.

    # 1. The kinetic term related to the imaginary time derivative
    kinetic_term_tau = sp.Derivative(phi, tau)**2
    
    # 2. The kinetic term related to the spatial gradient
    kinetic_term_space = grad_phi_sq
    
    # 3. The mass term (potential energy)
    mass_term = m**2 * phi**2
    
    # The complete integrand is the sum of these three terms.
    integrand = kinetic_term_tau + kinetic_term_space + mass_term
    
    # The inner product (ϕ, D_ϕ) is the integral of this expression over
    # (d+1)-dimensional Euclidean spacetime. The code will now print this integrand.
    
    print("The inner product (ϕ, D_ϕ) is the spacetime integral of the expression below.")
    print("The code has constructed this expression symbolically:")
    print("=" * 70)
    sp.pprint(integrand, use_unicode=True)
    print("=" * 70)
    print("This expression represents the integrand in the full inner product:")
    print("  (ϕ, D_ϕ) = ∫dτ∫dᵈx [ (∂ϕ/∂τ)² + (∇ϕ)² + m²ϕ² ]")
    print("\nWhere:")
    print("  - ϕ: The neutral scalar field, a function of τ and spatial coordinates.")
    print("  - m: The mass of the field quanta. The '2' is an exponent.")
    print("  - (∇ϕ)²: The square of the spatial gradient of the field.")
    print("  - ∂ϕ/∂τ: The partial derivative of the field with respect to imaginary time.")

if __name__ == '__main__':
    display_inner_product_integrand()