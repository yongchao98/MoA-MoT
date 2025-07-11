import sympy

def solve_field_theory_inner_product():
    """
    This script uses the SymPy library to define and display the inner product
    (ϕ, D_ϕ) for a neutral scalar field in finite-temperature field theory.
    It shows the representation in both configuration space and momentum space.
    """
    # Initialize pretty printing for clear mathematical output
    sympy.init_printing(use_unicode=True)

    print("--- Inner Product in Configuration Space ---")

    # Define the necessary symbols for configuration space
    phi = sympy.Function('ϕ')
    tau = sympy.Symbol('τ', real=True) # Imaginary time
    x, y, z = sympy.symbols('x y z', real=True) # Spatial coordinates
    m = sympy.Symbol('m', real=True, positive=True) # Mass
    beta = sympy.Symbol('β', real=True, positive=True) # Inverse temperature, β = 1/T

    # Define the scalar field ϕ as a function of Euclidean spacetime
    phi_spacetime = phi(tau, x, y, z)

    # Define the Euclidean Klein-Gordon operator, D = -∂_τ² - ∇² + m²
    # The numbers in the operator definition are the exponents '2' and the implicit coefficients '-1' and '+1'.
    D_phi = -sympy.diff(phi_spacetime, tau, 2) \
            - sympy.diff(phi_spacetime, x, 2) \
            - sympy.diff(phi_spacetime, y, 2) \
            - sympy.diff(phi_spacetime, z, 2) \
            + m**2 * phi_spacetime

    print("The operator D, acting on the field ϕ, is D = -∂_τ² - ∇² + m².")
    print("The result Dϕ is:")
    sympy.pprint(D_phi)
    print("\nThe numbers appearing in the operator are the implicit coefficients -1 and 1, and the exponent 2.")

    # Define the inner product (ϕ, D_ϕ) as the integral of ϕ * (Dϕ)
    integrand = phi_spacetime * D_phi
    inner_product_integral = sympy.Integral(
        integrand,
        (tau, 0, beta),
        (x, -sympy.oo, sympy.oo),
        (y, -sympy.oo, sympy.oo),
        (z, -sympy.oo, sympy.oo)
    )

    print("\nThe full inner product (ϕ, D_ϕ) is the integral over Euclidean spacetime:")
    sympy.pprint(inner_product_integral)


    print("\n\n--- Inner Product in Momentum Space ---")
    # This representation is often more useful for calculations.

    # Define symbols for momentum space
    phi_tilde = sympy.Function('~ϕ') # Fourier mode of the field
    omega_n = sympy.Symbol('ω_n', real=True) # Matsubara frequency
    px, py, pz = sympy.symbols('p_x p_y p_z', real=True) # 3-momentum components
    p_squared = px**2 + py**2 + pz**2 # Squared magnitude of 3-momentum
    n = sympy.Symbol('n', integer=True) # Matsubara mode index
    pi = sympy.pi # The number pi

    # In momentum space, the differential operator D becomes an algebraic factor
    momentum_space_factor = (omega_n**2 + p_squared + m**2)

    print("In momentum space, the inner product becomes a sum over Matsubara modes (n) and an")
    print("integral over 3-momentum (p). The operator D becomes a multiplicative factor:")
    sympy.pprint(momentum_space_factor)

    # The Matsubara frequencies contain the numbers 2 and π
    matsubara_freq_def = (2 * pi * n) / beta

    print("\nThis factor includes the Matsubara frequency ω_n, which is defined as:")
    sympy.pprint(matsubara_freq_def)
    print(f"The numbers in this definition are 2 and π (approx. {float(pi):.4f}).")

    # The full term in the momentum space sum/integral
    momentum_term = momentum_space_factor * sympy.Abs(phi_tilde(omega_n, px, py, pz))**2
    print("\nThe term inside the sum and integral for (ϕ, D_ϕ) is:")
    sympy.pprint(momentum_term)
    print("\nThe complete expression is (1/β) * Σ_n ∫ [d³p / (2π)³] multiplied by the term above.")

if __name__ == '__main__':
    solve_field_theory_inner_product()