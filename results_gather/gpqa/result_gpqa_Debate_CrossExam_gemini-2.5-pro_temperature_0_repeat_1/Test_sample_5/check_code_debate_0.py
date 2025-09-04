import sympy as sp

def check_quantum_energy_spectrum():
    """
    This function verifies the derivation of the energy spectrum for the given potential.
    It uses symbolic mathematics to perform the coordinate transformation and energy calculation.
    """
    try:
        # Define symbolic variables used in the problem
        # Assume k, m, r are positive real numbers
        k, m, r = sp.symbols('k m r', positive=True)
        # theta is a real number
        theta = sp.symbols('theta', real=True)
        # x and y are Cartesian coordinates
        x, y = sp.symbols('x y')
        # hbar is Planck's constant (positive)
        hbar = sp.symbols('hbar', positive=True)
        # n_x and n_y are non-negative integer quantum numbers
        n_x, n_y = sp.symbols('n_x n_y', integer=True, nonneg=True)

        # --- Step 1: Verify the coordinate transformation of the potential ---
        # The potential given in polar coordinates
        V_polar = sp.Rational(1, 2) * k * r**2 + sp.Rational(3, 2) * k * r**2 * sp.cos(theta)**2

        # The transformation rules are x = r*cos(theta) and r^2 = x^2 + y^2
        # We can substitute r*cos(theta) with x and r^2 with x^2 + y^2
        # sympy is smart enough to see r**2 * cos(theta)**2 as (r*cos(theta))**2
        V_cartesian_derived = V_polar.subs({r**2: x**2 + y**2, r*sp.cos(theta): x})
        
        # The potential claimed in the solution's derivation
        V_cartesian_solution = 2 * k * x**2 + sp.Rational(1, 2) * k * y**2

        # Check if the derived Cartesian potential matches the one in the solution
        if sp.simplify(V_cartesian_derived - V_cartesian_solution) != 0:
            return (f"Incorrect coordinate transformation. "
                    f"The potential V(r,θ) transforms to V(x,y) = {sp.simplify(V_cartesian_derived)}, "
                    f"but the solution claims it is {V_cartesian_solution}.")

        # --- Step 2: Verify the effective spring constants and frequencies ---
        # The potential is a 2D anisotropic harmonic oscillator: V(x,y) = 1/2*k_x*x^2 + 1/2*k_y*y^2
        # From V(x,y) = 2*k*x^2 + 1/2*k*y^2, we identify k_x and k_y
        k_x = 4 * k  # because 1/2 * k_x * x^2 = 2*k*x^2
        k_y = k      # because 1/2 * k_y * y^2 = 1/2*k*y^2

        # The angular frequencies are omega = sqrt(k_eff / m)
        omega_x = sp.sqrt(k_x / m)
        omega_y = sp.sqrt(k_y / m)

        # --- Step 3: Verify the total energy calculation ---
        # The energy for a 1D quantum harmonic oscillator is E_n = (n + 1/2)*hbar*omega
        E_x = (n_x + sp.Rational(1, 2)) * hbar * omega_x
        E_y = (n_y + sp.Rational(1, 2)) * hbar * omega_y
        
        # The total energy is the sum of the individual energies
        E_total_calculated = sp.expand(E_x + E_y)

        # --- Step 4: Compare the calculated energy with the answer from option B ---
        # The expression from option B
        E_answer_B = (2*n_x + n_y + sp.Rational(3, 2)) * hbar * sp.sqrt(k/m)
        E_answer_B_expanded = sp.expand(E_answer_B)

        # Check if the calculated total energy matches the expression from option B
        if sp.simplify(E_total_calculated - E_answer_B_expanded) != 0:
            return (f"Final energy expression is incorrect. "
                    f"The calculated energy is E = {E_total_calculated}, "
                    f"but the answer B provides E = {E_answer_B_expanded}.")

        # If all checks pass, the derivation and the final answer are correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_quantum_energy_spectrum()
print(result)