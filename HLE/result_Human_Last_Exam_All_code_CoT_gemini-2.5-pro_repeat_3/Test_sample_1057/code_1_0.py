import sympy

def solve_joule_heating_sphere():
    """
    Calculates the Joule heat dissipated by a shrinking, charge-leaking sphere.

    The solution is derived symbolically using sympy.
    """
    # Define the physical constants and variables as symbolic objects
    a = sympy.Symbol('a', positive=True, real=True)      # Initial radius
    V = sympy.Symbol('V', real=True)                    # Initial potential
    epsilon_0 = sympy.Symbol('epsilon_0', positive=True, real=True) # Vacuum permittivity
    pi = sympy.pi

    # --- Step 1: Calculate the initial electrostatic energy (U_initial) ---
    # The potential of a sphere is V = Q / (4*pi*epsilon_0*r).
    # Therefore, the initial charge Q_i on the sphere of radius 'a' at potential 'V' is:
    Q_initial = 4 * pi * epsilon_0 * a * V

    # The electrostatic energy of a charged sphere is U = (1/2) * Q * V.
    U_initial = sympy.Rational(1, 2) * Q_initial * V

    # --- Step 2: Analyze the energy conversion ---
    # The initial energy U_initial is converted into Joule heat (Q_Joule) and
    # mechanical work (W_mech) done by the electrostatic force.
    # U_initial = Q_Joule + W_mech
    
    # --- Step 3: Apply the "slowly shrinking" interpretation ---
    # The problem states the radius decreases "slowly". We interpret this to mean
    # the charge leakage happens much faster than the shrinking process.
    # Thus, we can assume the charge leaks completely while the radius is still 'a'.
    # In this scenario, the change in radius during leakage is zero (dr=0).
    # Mechanical work W_mech = integral(F_elec * dr) is therefore zero.
    
    W_mech = 0
    Q_Joule = U_initial - W_mech
    
    # --- Step 4: Print the derivation and the final result ---
    print("Derivation of the Joule Heat Dissipated:")
    print("------------------------------------------")
    print(f"1. The initial charge on the sphere is Q_i = 4 \u00b7 \u03c0 \u00b7 \u03b5\u2080 \u00b7 a \u00b7 V.")
    print(f"2. The initial stored electrostatic energy is U_i = (1/2) \u00b7 Q_i \u00b7 V.")
    print(f"   Substituting Q_i, we get U_i = {sympy.pretty(U_initial, use_unicode=True)}")
    
    print("\n3. This energy is converted into Joule heat (Q_Joule) and mechanical work (W_mech).")
    print("   The problem states the radius decreases 'slowly', implying that charge leakage is a much faster process.")
    print("   Therefore, we can assume the charge leaks completely while the radius is effectively constant at 'a'.")
    print("   Under this condition, the mechanical work done by the field is zero (W_mech = 0).")
    
    print("\n4. Thus, all the initial electrostatic energy is dissipated as Joule heat.")
    
    print("\nFinal equation for the Joule heat:")
    # Using unicode characters for a clean mathematical representation
    final_equation_unicode = f"Q_Joule = {sympy.pretty(Q_Joule, use_unicode=True)}"
    print(final_equation_unicode)

    # As requested, output the numbers in the final equation.
    # The equation is Q_Joule = 2 * pi * epsilon_0 * a * V**2
    print("\nNumerical components of the final equation:")
    print("  - The coefficient is 2.")
    print("  - The exponent of the potential V is 2.")

if __name__ == '__main__':
    solve_joule_heating_sphere()