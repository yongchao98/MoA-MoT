import sympy

def check_quantum_energy_spectrum():
    """
    This function checks the correctness of the given answer for the energy spectrum
    of a 2D quantum particle in a specific potential.

    The steps are:
    1. Define the potential in polar coordinates.
    2. Convert the potential to Cartesian coordinates.
    3. Identify the potential as two uncoupled 1D quantum harmonic oscillators (QHOs).
    4. Determine the angular frequencies (ω_x, ω_y) for each oscillator.
    5. Calculate the total energy spectrum E = E_x + E_y.
    6. Compare the derived energy spectrum with the provided answer.
    """
    try:
        # Define symbolic variables
        r, theta, k, m, hbar = sympy.symbols('r theta k m hbar', positive=True, real=True)
        x, y = sympy.symbols('x y', real=True)
        n_x, n_y = sympy.symbols('n_x n_y', integer=True, nonneg=True)

        # Step 1 & 2: Define potential and convert to Cartesian coordinates
        # V(r, θ) = 1/2*k*r^2 + 3/2*k*r^2*cos^2(θ)
        # Using x = r*cos(θ) and r^2 = x^2 + y^2
        V_cartesian = sympy.simplify((1/2)*k*(x**2 + y**2) + (3/2)*k*x**2)
        
        # Expected V_cartesian = 2*k*x**2 + 1/2*k*y**2
        V_x_term = V_cartesian.coeff(x**2) * x**2
        V_y_term = V_cartesian.coeff(y**2) * y**2

        # Step 3 & 4: Determine angular frequencies from the standard QHO potential V = 1/2*m*ω^2*z^2
        # For x-direction: 1/2*m*ω_x^2 = 2*k
        omega_x_sq = sympy.solve(sympy.Eq( (1/2)*m*sympy.Symbol('omega_x')**2, V_x_term.coeff(x**2) ), sympy.Symbol('omega_x')**2)[0]
        omega_x = sympy.sqrt(omega_x_sq)

        # For y-direction: 1/2*m*ω_y^2 = 1/2*k
        omega_y_sq = sympy.solve(sympy.Eq( (1/2)*m*sympy.Symbol('omega_y')**2, V_y_term.coeff(y**2) ), sympy.Symbol('omega_y')**2)[0]
        omega_y = sympy.sqrt(omega_y_sq)

        # Step 5: Calculate the total energy spectrum
        # E = E_x + E_y = (n_x + 1/2)ħω_x + (n_y + 1/2)ħω_y
        E_derived = (n_x + sympy.S(1)/2) * hbar * omega_x + (n_y + sympy.S(1)/2) * hbar * omega_y
        E_derived_simplified = sympy.simplify(E_derived)

        # Step 6: Compare with the provided answer (B)
        # Answer B: E = (2*n_x + n_y + 3/2) * ħ * sqrt(k/m)
        E_answer_B = (2*n_x + n_y + sympy.S(3)/2) * hbar * sympy.sqrt(k/m)

        # Check if the derived expression is equivalent to the answer's expression
        if sympy.simplify(E_derived_simplified - E_answer_B) == 0:
            return "Correct"
        else:
            return (f"Incorrect. The derived energy spectrum is E = {E_derived_simplified}, "
                    f"but the provided answer is E = {E_answer_B}. The derivation shows that "
                    f"the potential in Cartesian coordinates is V(x,y) = {V_cartesian}, leading to "
                    f"angular frequencies ω_x = {omega_x} and ω_y = {omega_y}. The sum of the "
                    f"1D QHO energies E_x + E_y results in the derived expression, which does not match the answer.")

    except Exception as e:
        return f"An error occurred during the checking process: {e}"

# Execute the check and print the result
result = check_quantum_energy_spectrum()
print(result)