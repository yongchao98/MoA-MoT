import sympy

def check_quantum_energy_spectrum():
    """
    Symbolically derives the energy spectrum for the given potential and checks
    the correctness of the provided answer.
    """
    try:
        # Define symbolic variables used in the problem
        m, k, hbar = sympy.symbols('m k hbar', positive=True)
        r, theta = sympy.symbols('r theta')
        x, y = sympy.symbols('x y')
        n_x, n_y = sympy.symbols('n_x n_y', integer=True, nonneg=True)

        # --- Step 1: Convert Potential to Cartesian Coordinates ---
        # The potential is given in polar coordinates
        V_polar = sympy.Rational(1, 2) * k * r**2 + sympy.Rational(3, 2) * k * r**2 * sympy.cos(theta)**2

        # Perform the substitution using r^2 = x^2 + y^2 and x = r*cos(theta)
        V_cartesian = V_polar.subs({
            r**2: x**2 + y**2,
            r*sympy.cos(theta): x
        })
        
        # Simplify the resulting expression
        V_cartesian_simplified = sympy.simplify(V_cartesian)
        
        # The expected potential in Cartesian coordinates is V(x,y) = 2kx^2 + 1/2 ky^2
        V_expected = 2 * k * x**2 + sympy.Rational(1, 2) * k * y**2
        
        if V_cartesian_simplified != V_expected:
            return f"Constraint Failure: Potential conversion to Cartesian coordinates is incorrect. Expected {V_expected}, but derived {V_cartesian_simplified}."

        # --- Step 2 & 3: Identify Components and Calculate Frequencies ---
        # The standard 1D Quantum Harmonic Oscillator (QHO) potential is V(z) = 1/2 * m * omega^2 * z^2
        
        # For the x-direction: V_x(x) = 2kx^2
        # Equating: 1/2 * m * omega_x^2 = 2k  => omega_x = 2*sqrt(k/m)
        V_x_coeff = V_cartesian_simplified.coeff(x**2)
        omega_x_sq_derived = 2 * V_x_coeff / m
        omega_x_derived = sympy.sqrt(omega_x_sq_derived)
        omega_x_expected = 2 * sympy.sqrt(k/m)
        
        if sympy.simplify(omega_x_derived - omega_x_expected) != 0:
            return f"Constraint Failure: Calculation of angular frequency ω_x is incorrect. Expected {omega_x_expected}, but derived {omega_x_derived}."

        # For the y-direction: V_y(y) = 1/2 ky^2
        # Equating: 1/2 * m * omega_y^2 = 1/2 k => omega_y = sqrt(k/m)
        V_y_coeff = V_cartesian_simplified.coeff(y**2)
        omega_y_sq_derived = 2 * V_y_coeff / m
        omega_y_derived = sympy.sqrt(omega_y_sq_derived)
        omega_y_expected = sympy.sqrt(k/m)

        if sympy.simplify(omega_y_derived - omega_y_expected) != 0:
            return f"Constraint Failure: Calculation of angular frequency ω_y is incorrect. Expected {omega_y_expected}, but derived {omega_y_derived}."

        # --- Step 4: Calculate Total Energy Spectrum ---
        # The energy for a 1D QHO is E_n = (n + 1/2) * hbar * omega
        E_x = (n_x + sympy.Rational(1, 2)) * hbar * omega_x_derived
        E_y = (n_y + sympy.Rational(1, 2)) * hbar * omega_y_derived
        
        # The total energy is the sum of the two
        E_total_derived = sympy.simplify(E_x + E_y)

        # --- Step 5: Compare with the Given Answer ---
        # The provided answer is A, which corresponds to the formula:
        # E = (2*n_x + n_y + 3/2) * hbar * sqrt(k/m)
        final_answer_formula = (2*n_x + n_y + sympy.Rational(3, 2)) * hbar * sympy.sqrt(k/m)
        
        # Check if the derived total energy matches the formula from the answer
        if sympy.simplify(E_total_derived - final_answer_formula) == 0:
            return "Correct"
        else:
            return (f"Incorrect: The final formula does not match the derivation.\n"
                    f"Derived Energy: {E_total_derived}\n"
                    f"Answer's Formula (Option A): {final_answer_formula}")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check
result = check_quantum_energy_spectrum()
print(result)