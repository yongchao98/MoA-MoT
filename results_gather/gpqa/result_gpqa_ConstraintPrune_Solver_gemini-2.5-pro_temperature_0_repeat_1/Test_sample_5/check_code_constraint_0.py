import sympy

def check_quantum_energy_spectrum():
    """
    This function verifies the energy spectrum for the given 2D quantum potential.
    It uses symbolic mathematics to perform the derivation and checks if the result
    matches the provided answer (Option D).
    """
    # 1. Define all necessary symbolic variables
    # Physical constants and parameters
    k, m, hbar = sympy.symbols('k m hbar', real=True, positive=True)
    # Polar coordinates
    r, theta = sympy.symbols('r theta', real=True)
    # Cartesian coordinates
    x, y = sympy.symbols('x y', real=True)
    # Quantum numbers
    n_x, n_y = sympy.symbols('n_x n_y', integer=True, nonneg=True)

    # 2. Define the potential in polar coordinates as given in the question
    V_polar = sympy.Rational(1, 2) * k * r**2 + sympy.Rational(3, 2) * k * r**2 * sympy.cos(theta)**2

    # 3. Convert the potential to Cartesian coordinates
    # Transformation rules: r^2 = x^2 + y^2 and x = r*cos(theta)
    # We can substitute r^2*cos(theta)^2 with x^2 and r^2 with x^2+y^2
    V_cartesian = V_polar.subs({
        r**2 * sympy.cos(theta)**2: x**2,
        r**2: x**2 + y**2
    })
    V_cartesian = sympy.simplify(V_cartesian)

    # Expected form: V(x,y) = 2*k*x**2 + 1/2*k*y**2
    expected_V_cartesian = 2 * k * x**2 + sympy.Rational(1, 2) * k * y**2
    if sympy.simplify(V_cartesian - expected_V_cartesian) != 0:
        return f"Constraint failed: Potential conversion to Cartesian coordinates is incorrect. Derived V(x,y) = {V_cartesian}."

    # 4. Identify the harmonic oscillator frequencies (omega_x, omega_y)
    # The general form is V(q) = 1/2 * m * omega^2 * q^2
    # For the x-direction: 2*k*x^2 = 1/2 * m * omega_x^2 * x^2  => omega_x^2 = 4*k/m
    # For the y-direction: 1/2*k*y^2 = 1/2 * m * omega_y^2 * y^2 => omega_y^2 = k/m
    omega_x = sympy.sqrt(4 * k / m)
    omega_y = sympy.sqrt(k / m)

    # 5. Calculate the total energy spectrum
    # The energy for a 1D QHO is E_n = (n + 1/2)*hbar*omega
    # The total energy is the sum of the energies for each independent dimension.
    E_derived = (n_x + sympy.Rational(1, 2)) * hbar * omega_x + (n_y + sympy.Rational(1, 2)) * hbar * omega_y
    E_derived_simplified = sympy.simplify(E_derived)

    # 6. Define the expression from the selected answer (Option D)
    E_option_D = (2 * n_x + n_y + sympy.Rational(3, 2)) * hbar * sympy.sqrt(k / m)

    # 7. Compare the derived energy with the expression from Option D
    # The difference should simplify to zero if they are identical.
    if sympy.simplify(E_derived_simplified - E_option_D) == 0:
        return "Correct"
    else:
        return (f"Incorrect: The derived energy spectrum is E = {E_derived_simplified}, "
                f"which does not match the expression from option D: E = {E_option_D}.")

# Execute the check and print the result
result = check_quantum_energy_spectrum()
print(result)