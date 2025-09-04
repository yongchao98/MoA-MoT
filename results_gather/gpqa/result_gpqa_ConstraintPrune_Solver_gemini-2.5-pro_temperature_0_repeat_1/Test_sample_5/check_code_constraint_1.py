import sympy

def check_quantum_energy_spectrum():
    """
    This function verifies the solution to the quantum mechanics problem by:
    1. Converting the potential from polar to Cartesian coordinates.
    2. Identifying the system as a 2D anisotropic harmonic oscillator.
    3. Deriving the angular frequencies for the x and y dimensions.
    4. Constructing the total energy spectrum from the derived frequencies.
    5. Comparing the derived spectrum with the provided answer (Option D).
    """
    try:
        # Define the necessary symbolic variables
        k, m, r, theta, x, y, hbar = sympy.symbols('k m r theta x y hbar', real=True, positive=True)
        n_x, n_y = sympy.symbols('n_x n_y', integer=True, nonneg=True)

        # --- Step 1: Convert potential to Cartesian coordinates ---
        # The potential is V(r, θ) = 1/2*k*r^2 + 3/2*k*r^2*cos^2(θ)
        # We can write this as V = 1/2*k*r^2 + 3/2*k*(r*cos(θ))^2
        # Using substitutions r^2 = x^2 + y^2 and r*cos(θ) = x
        V_cartesian = sympy.expand(
            (sympy.Rational(1, 2) * k * (x**2 + y**2)) + (sympy.Rational(3, 2) * k * x**2)
        )
        
        # Expected form: V(x,y) = 2*k*x**2 + 1/2*k*y**2
        expected_V = 2*k*x**2 + sympy.Rational(1, 2)*k*y**2
        if sympy.simplify(V_cartesian - expected_V) != 0:
            return f"Constraint Failure: Potential conversion to Cartesian coordinates is incorrect. Got {V_cartesian}, expected {expected_V}."

        # --- Step 2: Find the angular frequencies ω_x and ω_y ---
        # The general form of a 2D harmonic oscillator potential is V = 1/2*m*ω_x^2*x^2 + 1/2*m*ω_y^2*y^2
        # By comparing coefficients with our V_cartesian:
        # For x: 1/2*m*ω_x^2 = 2*k  => ω_x = 2*sqrt(k/m)
        # For y: 1/2*m*ω_y^2 = 1/2*k => ω_y = sqrt(k/m)
        
        omega_x = 2 * sympy.sqrt(k/m)
        omega_y = sympy.sqrt(k/m)

        # --- Step 3: Construct the total energy spectrum ---
        # The energy spectrum for a 2D anisotropic QHO is E = (n_x + 1/2)ħω_x + (n_y + 1/2)ħω_y
        derived_energy = (n_x + sympy.Rational(1, 2)) * hbar * omega_x + \
                         (n_y + sympy.Rational(1, 2)) * hbar * omega_y
        
        # Simplify the derived energy expression
        derived_energy_simplified = sympy.simplify(derived_energy)

        # --- Step 4: Compare with the provided answer (Option D) ---
        # Option D is E = (2*n_x + n_y + 3/2)ħ*sqrt(k/m)
        answer_d_expr = (2*n_x + n_y + sympy.Rational(3, 2)) * hbar * sympy.sqrt(k/m)

        # Check if the derived energy matches the expression from option D
        if sympy.simplify(derived_energy_simplified - answer_d_expr) == 0:
            return "Correct"
        else:
            return (f"Incorrect: The derived energy spectrum does not match the provided answer D.\n"
                    f"Derived Energy: {derived_energy_simplified}\n"
                    f"Answer D's Expression: {answer_d_expr}")

    except Exception as e:
        return f"An error occurred during verification: {e}"

# Execute the check
result = check_quantum_energy_spectrum()
print(result)