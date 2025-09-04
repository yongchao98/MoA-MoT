import sympy as sp

def check_quantum_oscillator_energy():
    """
    Verifies the energy spectrum for a 2D quantum harmonic oscillator with a given potential.

    The function follows these steps:
    1. Defines the potential in polar coordinates and converts it to Cartesian coordinates.
    2. Identifies the potential as a separable 2D anisotropic harmonic oscillator.
    3. Derives the angular frequencies (omega_x, omega_y) for the x and y directions.
    4. Constructs the total energy spectrum E = E_x + E_y.
    5. Compares the derived energy spectrum with the provided answer (Option D).
    """
    try:
        # Define symbolic variables used in the problem
        m, k, hbar = sp.symbols('m k hbar', positive=True, real=True)
        r, theta = sp.symbols('r theta', real=True)
        x, y = sp.symbols('x y', real=True)
        n_x, n_y = sp.symbols('n_x n_y', integer=True, nonnegative=True)

        # --- Step 1: Convert potential from polar to Cartesian coordinates ---
        # The potential is V(r, θ) = 1/2*k*r^2 + 3/2*k*r^2*cos^2(θ)
        # We can rewrite this as V = 1/2*k*r^2 + 3/2*k*(r*cos(θ))^2
        # Using x = r*cos(θ) and r^2 = x^2 + y^2
        V_cartesian_derived = sp.simplify(sp.Rational(1, 2) * k * (x**2 + y**2) + sp.Rational(3, 2) * k * x**2)
        
        # The expected form of the potential in Cartesian coordinates
        V_cartesian_expected = 2 * k * x**2 + sp.Rational(1, 2) * k * y**2
        
        if sp.simplify(V_cartesian_derived - V_cartesian_expected) != 0:
            return f"Incorrect potential conversion. Derived V(x,y) = {V_cartesian_derived}, but it should be {V_cartesian_expected}."

        # --- Step 2: Identify oscillator frequencies ---
        # The general form is V(x,y) = 1/2*m*omega_x^2*x^2 + 1/2*m*omega_y^2*y^2
        # Comparing coefficients for the x-term: 1/2*m*omega_x^2 = 2*k
        omega_x_sq = (2 * k) / (sp.Rational(1, 2) * m)
        omega_x = sp.sqrt(omega_x_sq)

        # Comparing coefficients for the y-term: 1/2*m*omega_y^2 = 1/2*k
        omega_y_sq = (sp.Rational(1, 2) * k) / (sp.Rational(1, 2) * m)
        omega_y = sp.sqrt(omega_y_sq)

        # Check if frequencies are derived correctly
        if sp.simplify(omega_x - 2*sp.sqrt(k/m)) != 0:
            return f"Incorrect frequency for x-direction. Derived omega_x = {omega_x}, but it should be 2*sqrt(k/m)."
        if sp.simplify(omega_y - sp.sqrt(k/m)) != 0:
            return f"Incorrect frequency for y-direction. Derived omega_y = {omega_y}, but it should be sqrt(k/m)."

        # --- Step 3: Construct the total energy spectrum ---
        # E = E_x + E_y = (n_x + 1/2)ħω_x + (n_y + 1/2)ħω_y
        E_derived = (n_x + sp.Rational(1, 2)) * hbar * omega_x + (n_y + sp.Rational(1, 2)) * hbar * omega_y
        E_derived_simplified = sp.simplify(E_derived)

        # --- Step 4: Compare with the given answer (Option D) ---
        # Answer D is E = (2*n_x + n_y + 3/2)ħ*sqrt(k/m)
        E_answer_D = (2*n_x + n_y + sp.Rational(3, 2)) * hbar * sp.sqrt(k/m)

        # The difference between the derived energy and the answer should be zero
        if sp.simplify(E_derived_simplified - E_answer_D) != 0:
            return (f"The derived energy spectrum does not match the answer D.\n"
                    f"Derived Energy: {E_derived_simplified}\n"
                    f"Answer D Energy: {E_answer_D}")

        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_quantum_oscillator_energy()
print(result)