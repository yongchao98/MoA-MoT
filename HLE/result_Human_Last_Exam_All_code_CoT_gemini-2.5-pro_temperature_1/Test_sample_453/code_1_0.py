import sympy

def calculate_energy_shift():
    """
    Calculates the relativistic energy shift for a hydrogen atom state symbolically.
    """
    # Define fundamental constants and quantum numbers as symbols
    m, e, hbar, c, epsilon_0 = sympy.symbols('m e ħ c ε₀', real=True, positive=True)
    n, ell = sympy.symbols('n ℓ', integer=True, positive=True)

    # Given quantum numbers for the state
    n_val = 3
    ell_val = 2

    # --- Step 1: Define key expressions ---
    # Coulomb constant k_e in terms of fundamental constants
    k_e = 1 / (4 * sympy.pi * epsilon_0)
    # Bohr radius a_0
    a_0 = (4 * sympy.pi * epsilon_0 * hbar**2) / (m * e**2)
    # Unperturbed energy E_n^(0)
    E_n0 = - (k_e * e**2) / (2 * a_0 * n**2)

    # --- Step 2: Calculate expectation values ---
    # <V> from Virial Theorem: <V> = 2 * E_n^(0)
    V_exp = 2 * E_n0
    # <1/r^2> formula for hydrogen atom
    inv_r2_exp = 1 / (n**3 * (ell + sympy.S(1)/2) * a_0**2)
    # <V^2> = (k_e * e^2)^2 * <1/r^2>
    V2_exp = (k_e * e**2)**2 * inv_r2_exp

    # --- Step 3: Calculate <p^4> ---
    # <p^4> = 4*m^2 * ( (E_n0)^2 - 2*E_n0*<V> + <V^2> )
    p4_exp = 4 * m**2 * sympy.simplify(E_n0**2 - 2 * E_n0 * V_exp + V2_exp)

    # --- Step 4: Calculate the final energy shift ΔE ---
    # ΔE = -<p^4> / (8 * m^3 * c^2)
    Delta_E = -p4_exp / (8 * m**3 * c**2)

    # --- Step 5: Substitute n=3, l=2 and simplify ---
    final_expression = Delta_E.subs({n: n_val, ell: ell_val})
    simplified_expression = sympy.simplify(final_expression)

    # --- Step 6: Print the result in a readable format ---
    # The expected form is -C * m * e^8 / ( (4*pi*epsilon_0)^4 * hbar^4 * c^2 )
    # Let's extract the coefficient C to present the equation clearly.
    num, den = simplified_expression.as_numer_denom()
    num_coeff = abs(num.as_coeff_mul()[0])
    den_coeff_part = den / ((4 * sympy.pi * epsilon_0)**4 * hbar**4 * c**2)
    den_coeff = sympy.simplify(den_coeff_part)

    print("The calculated energy shift due to the relativistic correction is:")
    # Using unicode characters for a clear representation
    print(f"\n  ΔE = - ( {num_coeff} ⋅ m ⋅ e⁸ ) / ( {den_coeff} ⋅ (4πε₀)⁴ ⋅ ħ⁴ ⋅ c² )")

    print("\nThe numerical values in this final equation are:")
    print(f"  - Coefficient in the numerator: {num_coeff}")
    print(f"  - Exponent of electron charge 'e': 8")
    print(f"  - Coefficient in the denominator: {den_coeff}")
    print(f"  - Exponent of the Coulomb term '(4πε₀)': 4")
    print(f"  - Exponent of reduced Planck constant 'ħ': 4")
    print(f"  - Exponent of speed of light 'c': 2")

if __name__ == '__main__':
    calculate_energy_shift()