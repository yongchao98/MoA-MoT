import math

def solve_joule_heat():
    """
    This function explains the derivation and prints the final formula for the
    Joule heat dissipated by the leaking and shrinking charged sphere.
    """
    # Symbolic variables (represented as strings for the explanation)
    a = "a"
    V = "V"
    epsilon_0 = "ε₀"
    pi = "π"

    print("Derivation of the Total Joule Heat Dissipated:")
    print("=" * 50)

    # Step 1: State the governing principle (Energy Conservation)
    print("1. The total Joule heat (H) dissipated is equal to the initial electrostatic potential energy (U_initial) stored in the sphere.")
    print("   This is because the final charge is zero, making the final energy zero.\n")

    # Step 2: Provide the formula for the energy of a charged sphere
    print(f"2. The initial energy is given by the formula: U_initial = (1/2) * C * {V}²")
    print(f"   where C is the initial capacitance and {V} is the initial potential.\n")

    # Step 3: Provide the formula for the capacitance of a sphere
    print(f"3. The capacitance (C) of a sphere with radius '{a}' is: C = 4 * {pi} * {epsilon_0} * {a}\n")

    # Step 4: Combine the formulas to find the initial energy
    print("4. Substituting the capacitance into the energy formula:")
    print(f"   U_initial = (1/2) * (4 * {pi} * {epsilon_0} * {a}) * {V}²")
    print("   Simplifying this gives the total Joule heat (H).\n")

    # Step 5: Present the final equation
    print("Final Equation for Joule Heat (H):")
    print("-" * 50)
    # The final equation explicitly shows the numbers involved.
    final_equation = f"H = 2 * {pi} * {epsilon_0} * {a} * {V}²"
    print(final_equation)
    print("\nComponent breakdown of the equation:")
    print(f"  - The number in the equation is: 2")
    print(f"  - {pi} is the mathematical constant Pi (approx. {math.pi:.5f})")
    print(f"  - {epsilon_0} is the permittivity of free space (a physical constant)")
    print(f"  - {a} is the initial radius of the sphere")
    print(f"  - {V} is the initial potential of the sphere")

# Execute the function to print the solution
solve_joule_heat()