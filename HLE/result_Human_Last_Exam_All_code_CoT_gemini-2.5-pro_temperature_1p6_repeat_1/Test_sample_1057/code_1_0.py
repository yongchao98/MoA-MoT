def solve_sphere_heat_dissipation():
    """
    This script derives and prints the formula for the total Joule heat
    dissipated by a charged sphere that leaks charge while shrinking.
    """
    # Using string representations for physical constants and variables
    # as no numerical values are provided.
    V = "V"
    a = "a"
    pi = "π"
    epsilon_0 = "ε₀"

    print("Derivation of the total Joule heat dissipated:")
    
    print("\nStep 1: The electrostatic potential energy (U) stored in the sphere is determined.")
    print("The formula for stored energy in a capacitor is U = (1/2) * C * V^2.")
    print(f"The capacitance (C) of a sphere with radius 'a' is C = 4 * {pi} * {epsilon_0} * {a}.")
    
    print("\nStep 2: The initial energy (U_initial) is calculated by substituting C.")
    print(f"U_initial = (1/2) * (4 * {pi} * {epsilon_0} * {a}) * {V}^2")

    print("\nStep 3: Simplifying the expression gives the initial stored energy.")
    # The coefficients are 1/2 and 4. The product is 2.
    factor = 2
    print(f"U_initial = {factor} * {pi} * {epsilon_0} * {a} * {V}^2")

    print("\nStep 4: The final Joule heat (H) is equal to the initial stored energy.")
    print("As the sphere loses all its charge, this initial energy is entirely converted into heat and mechanical work.")
    print("Assuming all the initial potential energy is dissipated as Joule heat, we get the final result.")

    print("\n-------------------------------------------")
    print("The total Joule heat (H) dissipated is:")
    print("H = ", end="")
    print(f"{factor} * ", end="")
    print(f"{pi} * ", end="")
    print(f"{epsilon_0} * ", end="")
    print(f"{a} * ", end="")
    print(f"{V}^2")
    print("-------------------------------------------")

solve_sphere_heat_dissipation()
