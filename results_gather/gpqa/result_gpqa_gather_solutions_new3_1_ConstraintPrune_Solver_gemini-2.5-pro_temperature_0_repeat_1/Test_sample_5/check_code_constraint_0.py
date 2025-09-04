import sympy

def check_correctness():
    """
    Checks the correctness of the provided answer for the energy spectrum.

    The verification process follows these steps:
    1.  Defines the potential in polar coordinates using symbolic variables.
    2.  Converts the potential to Cartesian coordinates (x, y).
    3.  Simplifies the Cartesian potential and verifies it matches the expected form of a 2D anisotropic harmonic oscillator.
    4.  Calculates the angular frequencies (omega_x, omega_y) based on the standard QHO potential form V = 1/2 * m * omega^2 * z^2.
    5.  Derives the total energy spectrum by summing the energies of the two independent 1D oscillators.
    6.  Compares the derived energy spectrum with the formula from the provided answer (Option D).
    """
    try:
        # Define symbolic variables
        # 'real=True, positive=True' for physical constants
        r, theta, k, m, hbar = sympy.symbols('r theta k m hbar', real=True, positive=True)
        x, y = sympy.symbols('x y', real=True)
        # 'integer=True, nonneg=True' for quantum numbers
        n_x, n_y = sympy.symbols('n_x n_y', integer=True, nonneg=True)

        # --- Step 1: Coordinate Transformation ---
        # Original potential in polar coordinates
        V_polar = sympy.S(1)/2 * k * r**2 + sympy.S(3)/2 * k * r**2 * sympy.cos(theta)**2

        # Transformation rules to Cartesian coordinates
        transformations = {
            r: sympy.sqrt(x**2 + y**2),
            sympy.cos(theta): x / sympy.sqrt(x**2 + y**2)
        }

        # Substitute to get potential in Cartesian coordinates and simplify
        V_cartesian = V_polar.subs(transformations)
        V_cartesian_simplified = sympy.simplify(V_cartesian)

        # Expected Cartesian potential from the derivation
        expected_V_cartesian = 2*k*x**2 + sympy.S(1)/2*k*y**2

        if sympy.simplify(V_cartesian_simplified - expected_V_cartesian) != 0:
            return f"Constraint Failure: The potential transformation to Cartesian coordinates is incorrect. Expected {expected_V_cartesian}, but the code derived {V_cartesian_simplified}."

        # --- Step 2: Identify Frequencies ---
        # The standard 1D QHO potential is V(z) = 1/2 * m * omega^2 * z^2
        # From V_cartesian, we identify the potential terms for x and y
        V_x_coeff = V_cartesian_simplified.coeff(x**2) # This should be 2*k
        V_y_coeff = V_cartesian_simplified.coeff(y**2) # This should be k/2

        # Define omega symbols
        omega_x, omega_y = sympy.symbols('omega_x omega_y', real=True, positive=True)

        # Solve for omega_x
        eq_x = sympy.Eq(sympy.S(1)/2 * m * omega_x**2, V_x_coeff)
        sol_omega_x = sympy.solve(eq_x, omega_x)[0]
        expected_omega_x = 2 * sympy.sqrt(k/m)
        if sympy.simplify(sol_omega_x - expected_omega_x) != 0:
            return f"Constraint Failure: The angular frequency for the x-direction is incorrect. Expected {expected_omega_x}, but got {sol_omega_x}."

        # Solve for omega_y
        eq_y = sympy.Eq(sympy.S(1)/2 * m * omega_y**2, V_y_coeff)
        sol_omega_y = sympy.solve(eq_y, omega_y)[0]
        expected_omega_y = sympy.sqrt(k/m)
        if sympy.simplify(sol_omega_y - expected_omega_y) != 0:
            return f"Constraint Failure: The angular frequency for the y-direction is incorrect. Expected {expected_omega_y}, but got {sol_omega_y}."

        # --- Step 3: Calculate Total Energy ---
        # Energy for 1D QHO: E_n = (n + 1/2) * hbar * omega
        E_x = (n_x + sympy.S(1)/2) * hbar * sol_omega_x
        E_y = (n_y + sympy.S(1)/2) * hbar * sol_omega_y
        E_total_derived = sympy.simplify(E_x + E_y)

        # --- Step 4: Compare with the Provided Answer ---
        # The provided answer's reasoning and final choice is D: E = (2*n_x + n_y + 3/2) * hbar * sqrt(k/m)
        answer_expression = (2*n_x + n_y + sympy.S(3)/2) * hbar * sympy.sqrt(k/m)

        # Check if the derived expression matches the answer's expression
        if sympy.simplify(E_total_derived - answer_expression) != 0:
            return f"Constraint Failure: The final derived energy spectrum does not match the provided answer's formula.\nDerived: {E_total_derived}\nProvided Answer Formula: {answer_expression}"

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during verification: {e}"

# Run the check
result = check_correctness()
print(result)