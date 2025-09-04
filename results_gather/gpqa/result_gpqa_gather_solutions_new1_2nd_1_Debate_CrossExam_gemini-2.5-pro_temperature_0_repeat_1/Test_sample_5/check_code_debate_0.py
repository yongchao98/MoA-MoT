import sympy

def check_correctness_of_quantum_energy_spectrum():
    """
    This function programmatically verifies the derivation of the energy spectrum for the given potential.
    It uses the symbolic math library sympy to perform the coordinate transformation,
    derive the frequencies of the resulting harmonic oscillators, and calculate the total energy.
    Finally, it compares the derived result with the expression from the selected answer (Option C).
    """
    # Define symbols for the variables in the problem
    m, k, hbar = sympy.symbols('m k hbar', positive=True)
    r, theta = sympy.symbols('r theta')
    x, y = sympy.symbols('x y')
    n_x, n_y = sympy.symbols('n_x n_y', integer=True, nonneg=True)
    omega_x, omega_y = sympy.symbols('omega_x omega_y', positive=True)

    # Step 1: Define the potential in polar coordinates from the question
    V_polar = sympy.Rational(1, 2) * k * r**2 + sympy.Rational(3, 2) * k * r**2 * sympy.cos(theta)**2

    # Step 2: Convert the potential to Cartesian coordinates
    # Use the direct relations: x^2 = r^2*cos(theta)^2 and r^2 = x^2 + y^2
    V_cartesian = V_polar.subs({
        r**2 * sympy.cos(theta)**2: x**2,
        r**2: x**2 + y**2
    })
    V_cartesian = sympy.expand(V_cartesian)

    # The expected Cartesian potential after simplification is 2*k*x**2 + 1/2*k*y**2
    V_expected = 2 * k * x**2 + sympy.Rational(1, 2) * k * y**2
    if sympy.simplify(V_cartesian - V_expected) != 0:
        return f"Incorrect potential conversion. Derived V(x,y) = {V_cartesian}, but expected {V_expected}."

    # Step 3: Separate the potential and derive the angular frequencies
    V_x_component = V_cartesian.subs(y, 0)
    V_y_component = V_cartesian.subs(x, 0)

    # The standard 1D QHO potential is V = 1/2 * m * omega^2 * q^2
    # Solve for omega_x
    eq_x = sympy.Eq(V_x_component, sympy.Rational(1, 2) * m * omega_x**2 * x**2)
    sol_omega_x = sympy.solve(eq_x, omega_x)
    derived_omega_x = [sol for sol in sol_omega_x if sol.is_positive][0] # Select the positive solution

    # Solve for omega_y
    eq_y = sympy.Eq(V_y_component, sympy.Rational(1, 2) * m * omega_y**2 * y**2)
    sol_omega_y = sympy.solve(eq_y, omega_y)
    derived_omega_y = [sol for sol in sol_omega_y if sol.is_positive][0] # Select the positive solution

    # Step 4: Construct the total energy spectrum from the derived frequencies
    # The energy for a 1D QHO is E_n = (n + 1/2) * hbar * omega
    E_x = (n_x + sympy.Rational(1, 2)) * hbar * derived_omega_x
    E_y = (n_y + sympy.Rational(1, 2)) * hbar * derived_omega_y
    
    derived_E = sympy.expand(E_x + E_y)
    derived_E = sympy.factor(derived_E) # Factor to get a clean, comparable form

    # Step 5: Define the expression from the provided answer (Option C)
    # C) E = (2n_x+n_y+3/2)ℏ*sqrt(k/m)
    answer_C_expression = (2 * n_x + n_y + sympy.Rational(3, 2)) * hbar * sympy.sqrt(k/m)

    # Step 6: Compare the programmatically derived energy with the answer's expression
    # The difference should simplify to zero if they are equivalent.
    if sympy.simplify(derived_E - answer_C_expression) == 0:
        return "Correct"
    else:
        return (f"The final answer is incorrect. "
                f"The programmatically derived energy spectrum is E = {derived_E}. "
                f"The energy spectrum from answer C is E = {answer_C_expression}. "
                f"These expressions are not equivalent.")

# Run the check and print the result
result = check_correctness_of_quantum_energy_spectrum()
print(result)