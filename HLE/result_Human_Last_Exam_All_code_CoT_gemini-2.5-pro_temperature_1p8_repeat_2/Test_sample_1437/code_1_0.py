import sympy

def solve_lindhard_function():
    """
    Calculates the normalized value of the Lindhard polarization function
    at zero frequency and zero momentum transfer for a 3D electron gas.
    """
    # Step 1: Define symbolic variables for the physical quantities.
    # m: electron mass, hbar: reduced Planck constant, E_F: Fermi energy.
    m, hbar, E_F = sympy.symbols('m hbar E_F', positive=True)

    # Step 2: Express the electron density 'n' as a function of the Fermi energy 'E_F'.
    # For a 3D homogeneous electron gas at T=0, n is proportional to E_F^(3/2).
    # n = (1 / (3*pi**2)) * (2*m*E_F / hbar**2)**(3/2)
    n_expression = (1 / (3 * sympy.pi**2)) * (2 * m * E_F / hbar**2)**(sympy.S(3)/2)

    # Step 3: Calculate the density of states at the Fermi level, g(E_F).
    # At T=0, the density of states g(E_F) is the derivative of n with respect to E_F.
    g_EF = sympy.diff(n_expression, E_F)

    # Step 4: Use the standard result from linear response theory.
    # The static, long-wavelength Lindhard function Pi_0(q=0, w=0) is -g(E_F).
    Pi0_00 = -g_EF

    # Step 5: Determine the universal numerical value.
    # The value of Pi_0(0,0) depends on the system's parameters. To obtain a universal
    # numerical value, we normalize it by the density of states at the Fermi level, g(E_F).
    normalized_value = Pi0_00 / g_EF

    # Step 6: Print the explanation and the final result.
    print("The Lindhard polarization function at zero momentum and frequency, Pi_0(0,0), is equal to the negative of the density of states at the Fermi level, -g(E_F).")
    print("To obtain a universal numerical value, we evaluate the ratio of these two quantities.")
    print("\nFinal Equation:")
    
    # We construct the final equation with numbers as requested.
    # The ratio simplifies to -1 / 1.
    numerator = -1
    denominator = 1
    result = int(normalized_value)

    print(f"Value = Pi_0(0,0) / g(E_F)")
    print(f"Value = ({numerator} * g(E_F)) / ({denominator} * g(E_F))")
    print(f"Value = {result}")

if __name__ == "__main__":
    solve_lindhard_function()