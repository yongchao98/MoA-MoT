import sympy

def solve_sphere_heat_dissipation():
    """
    Calculates and explains the total Joule heat dissipated by a shrinking, leaking charged sphere.
    """
    # Define the symbols used in the equations
    a, V, epsilon_0, pi = sympy.symbols('a V ε₀ π')

    print("Step 1: Define the initial electrostatic energy of the sphere.")
    # Capacitance of a sphere of radius 'a'
    C = 4 * pi * epsilon_0 * a
    # Initial electrostatic energy U = (1/2) * C * V^2
    U_initial = sympy.Rational(1, 2) * C * V**2

    print(f"The capacitance of the sphere is C = {C}.")
    print(f"The initial stored energy is U_initial = (1/2) * C * V^2 = {U_initial}.\n")

    print("Step 2: Apply the principle of conservation of energy.")
    print("The initial stored energy is converted into two forms:")
    print("  1. Direct Joule heat (H) from charge leaking through the atmosphere.")
    print("  2. Mechanical work (W_gas) done on the gas as the sphere shrinks.")
    print("The total energy balance is: U_initial = H + W_gas.\n")

    print("Step 3: Determine the total heat dissipated into the atmosphere.")
    print("The work done on the gas (W_gas) increases its energy. This gas is released into the")
    print("atmosphere, where its energy thermalizes, also becoming heat.")
    print("Therefore, the total heat dissipated into the atmosphere is the sum of these two components: Total Heat = H + W_gas.\n")

    print("Step 4: Final Result.")
    print("From the energy balance, the total heat dissipated is equal to the initial energy.")
    
    # The final equation for the total heat dissipated
    total_heat = U_initial
    
    # The prompt requires outputting each number/symbol in the final equation.
    # The final equation is: Total Heat = 2 * π * ε₀ * a * V**2
    
    print("\n--- Final Equation ---")
    print(f"Total Heat = {total_heat}")
    print("----------------------")
    
    print("\nComponents of the final equation:")
    print("Numeric Coefficient: 2")
    print(f"Mathematical Constant: {pi}")
    print(f"Physical Constant: {epsilon_0} (Permittivity of free space)")
    print(f"Initial Condition: {a} (Initial radius of the sphere)")
    print(f"Initial Condition: {V}**2 (Square of the initial potential)")

solve_sphere_heat_dissipation()