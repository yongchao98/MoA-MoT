import sympy

def solve_energy_shift():
    """
    Calculates the ground state energy shift for two interacting
    quantum harmonic oscillators using second-order perturbation theory.
    """
    # Define the physical constants and variables as symbolic variables.
    # All are positive real numbers.
    e, m, omega0, R, hbar = sympy.symbols('e m omega_0 R hbar', real=True, positive=True)

    # 1. Define the perturbation Hamiltonian H'.
    # The leading term in the dipole-dipole interaction is H' = -C * x1 * x2.
    # In the specified units, the interaction constant C is 2*e^2/R^3.
    C = 2 * e**2 / R**3

    # 2. Calculate the squared matrix element |<1,1|H'|0,0>|^2.
    # The matrix element <1|x|0> for a single QHO is sqrt(hbar / (2*m*omega_0)).
    # So, |<1|x|0>|^2 is hbar / (2*m*omega_0).
    x01_squared = hbar / (2 * m * omega0)
    
    # The matrix element M = <1,1|H'|0,0> = -C * <1|x1|0> * <1|x2|0>.
    # The squared matrix element |M|^2 is C^2 * |<1|x1|0>|^2 * |<1|x2|0>|^2.
    M_squared = C**2 * x01_squared * x01_squared

    # 3. Define the energy denominator.
    # E_ground = hbar*omega_0.
    # E_excited (state |1,1>) = 3*hbar*omega_0.
    # The denominator is E_ground - E_excited.
    energy_denominator = hbar * omega0 - 3 * hbar * omega0

    # 4. Calculate the second-order energy shift.
    # dE_0 = |M|^2 / (E_ground - E_excited)
    delta_E0 = M_squared / energy_denominator

    # 5. Simplify the expression and print the final equation.
    final_expression = sympy.simplify(delta_E0)
    
    print("The final equation for the ground state energy shift (ΔE_0) is:")
    # The following print statement shows the final equation with all its numeric
    # coefficients and powers, as requested.
    print(f"ΔE_0 = {final_expression}")

solve_energy_shift()