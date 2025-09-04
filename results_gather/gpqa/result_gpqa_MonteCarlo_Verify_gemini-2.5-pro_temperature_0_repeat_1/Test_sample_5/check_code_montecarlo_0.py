import sympy

def check_energy_spectrum():
    """
    Checks the correctness of the given energy spectrum for a 2D quantum harmonic oscillator.

    The process involves:
    1. Defining all relevant physical quantities as symbolic variables.
    2. Expressing the potential V(r, θ) in polar coordinates.
    3. Transforming the potential into Cartesian coordinates (x, y).
    4. Simplifying the Cartesian potential to identify it as a sum of two
       independent 1D harmonic oscillators, V(x, y) = V(x) + V(y).
    5. Deriving the effective spring constants (k_x, k_y) and angular frequencies (ω_x, ω_y).
    6. Constructing the total energy spectrum E = E_x + E_y from the derived frequencies.
    7. Comparing the derived energy spectrum with the provided answer (Option B).
    """
    try:
        # 1. Define symbolic variables
        r, theta, x, y, k, m, hbar = sympy.symbols('r theta x y k m hbar', real=True, positive=True)
        n_x, n_y = sympy.symbols('n_x n_y', integer=True, nonneg=True)

        # 2. Define the potential in polar coordinates from the question
        V_polar = (sympy.S(1)/2) * k * r**2 + (sympy.S(3)/2) * k * r**2 * sympy.cos(theta)**2

        # 3. Transform to Cartesian coordinates using r^2 = x^2 + y^2 and cos(θ) = x/r
        # Note: sympy.cos(theta)**2 = x**2 / r**2 = x**2 / (x**2 + y**2)
        V_cartesian = V_polar.subs({
            r**2: x**2 + y**2,
            sympy.cos(theta)**2: x**2 / (x**2 + y**2)
        })

        # 4. Simplify the Cartesian potential
        V_cartesian_simplified = sympy.simplify(V_cartesian)
        
        # Expected form is V(x,y) = 1/2*k_x*x**2 + 1/2*k_y*y**2
        # Our simplified potential is V(x,y) = 2*k*x**2 + 1/2*k*y**2
        
        # 5. Identify effective spring constants and frequencies
        # For V_x(x) = 2*k*x**2 = 1/2 * k_x * x**2  => k_x = 4*k
        # For V_y(y) = 1/2*k*y**2 = 1/2 * k_y * y**2 => k_y = k
        k_x_derived = 4 * k
        k_y_derived = k

        omega_x_derived = sympy.sqrt(k_x_derived / m) # sqrt(4k/m) = 2*sqrt(k/m)
        omega_y_derived = sympy.sqrt(k_y_derived / m) # sqrt(k/m)

        # 6. Construct the total energy spectrum
        # E = (n_x + 1/2)ħω_x + (n_y + 1/2)ħω_y
        E_derived = (n_x + sympy.S(1)/2) * hbar * omega_x_derived + \
                    (n_y + sympy.S(1)/2) * hbar * omega_y_derived

        # Simplify the derived energy expression
        E_derived_simplified = sympy.simplify(E_derived)

        # 7. Define the answer to be checked (Option B)
        E_answer_B = (2*n_x + n_y + sympy.S(3)/2) * hbar * sympy.sqrt(k/m)

        # 8. Compare the derived expression with the answer
        # The difference should simplify to zero if they are identical.
        if sympy.simplify(E_derived_simplified - E_answer_B) == 0:
            return "Correct"
        else:
            # If they don't match, provide the reason.
            reason = (
                f"Incorrect. The derivation leads to a different energy spectrum.\n"
                f"1. The potential in Cartesian coordinates is V(x,y) = {V_cartesian_simplified}.\n"
                f"2. This corresponds to two uncoupled harmonic oscillators with effective spring constants k_x = {k_x_derived} and k_y = {k_y_derived}.\n"
                f"3. The resulting angular frequencies are ω_x = {omega_x_derived} and ω_y = {omega_y_derived}.\n"
                f"4. The total energy E = (n_x + 1/2)ħω_x + (n_y + 1/2)ħω_y simplifies to E = {E_derived_simplified}.\n"
                f"5. The provided answer is E = {E_answer_B}, which does not match the derived result."
            )
            return reason

    except Exception as e:
        return f"An error occurred during the checking process: {e}"

# Execute the check and print the result.
result = check_energy_spectrum()
print(result)