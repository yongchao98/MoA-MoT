import sympy

def check_energy_spectrum():
    """
    This function checks the correctness of the derived energy spectrum for a quantum particle
    in a given 2D potential.

    The potential is V(r, θ) = 1/2 kr^2 + 3/2 kr^2 cos^2(θ).
    The expected answer is E = (2*n_x + n_y + 3/2) * hbar * sqrt(k/m), which corresponds to option D.
    """
    try:
        # Define symbolic variables
        k, m, hbar = sympy.symbols('k m hbar', positive=True)
        r, theta = sympy.symbols('r theta')
        x, y = sympy.symbols('x y')
        n_x, n_y = sympy.symbols('n_x n_y', integer=True, nonneg=True)

        # --- Step 1: Convert the Potential to Cartesian Coordinates ---
        # V(r, θ) = 1/2 kr^2 + 3/2 kr^2 cos^2(θ)
        # Transformations: r^2 = x^2 + y^2 and x = r*cos(θ)
        V_cartesian = sympy.Rational(1, 2) * k * (x**2 + y**2) + sympy.Rational(3, 2) * k * x**2
        V_simplified = sympy.simplify(V_cartesian)
        
        # Expected Cartesian potential: 2*k*x**2 + 1/2*k*y**2
        expected_V = 2*k*x**2 + sympy.Rational(1, 2)*k*y**2
        if sympy.simplify(V_simplified - expected_V) != 0:
            return f"Incorrect potential conversion. Derived V(x,y) = {V_simplified}, but expected {expected_V}."

        # --- Step 2: Identify the System and Separate the Variables ---
        # The potential is a 2D anisotropic harmonic oscillator.
        # V(x,y) = V_x(x) + V_y(y) where V_z(z) = 1/2 * k_eff * z^2
        
        # Extract coefficients for x^2 and y^2
        coeff_x2 = V_simplified.coeff(x**2)
        coeff_y2 = V_simplified.coeff(y**2)

        # --- Step 3: Calculate the Effective Spring Constants and Angular Frequencies ---
        # For x-direction: V_x(x) = 2*k*x^2 = 1/2 * k_x * x^2 => k_x = 4k
        k_x = 2 * coeff_x2
        omega_x = sympy.sqrt(k_x / m)
        
        # For y-direction: V_y(y) = 1/2*k*y^2 = 1/2 * k_y * y^2 => k_y = k
        k_y = 2 * coeff_y2
        omega_y = sympy.sqrt(k_y / m)

        # --- Step 4: Calculate the Total Energy Spectrum ---
        # E_n = (n + 1/2) * hbar * omega
        E_x = (n_x + sympy.Rational(1, 2)) * hbar * omega_x
        E_y = (n_y + sympy.Rational(1, 2)) * hbar * omega_y
        
        E_derived = sympy.expand(E_x + E_y)

        # --- Step 5: Compare with the Given Answer ---
        # The provided answer corresponds to option D: E = (2*n_x + n_y + 3/2) * hbar * sqrt(k/m)
        E_answer = (2*n_x + n_y + sympy.Rational(3, 2)) * hbar * sympy.sqrt(k/m)
        E_answer_expanded = sympy.expand(E_answer)

        # Check if the derived expression matches the answer's expression
        if sympy.simplify(E_derived - E_answer_expanded) == 0:
            return "Correct"
        else:
            return (f"Incorrect final energy expression.\n"
                    f"Derived expression: {E_derived}\n"
                    f"Answer's expression: {E_answer_expanded}\n"
                    f"The derivation steps seem correct, but the final algebraic simplification or comparison failed.")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_energy_spectrum()
print(result)