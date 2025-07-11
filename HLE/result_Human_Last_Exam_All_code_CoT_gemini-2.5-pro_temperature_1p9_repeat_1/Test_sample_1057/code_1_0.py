def solve_sphere_heat_dissipation():
    """
    This function derives and prints the formula for the total Joule heat
    dissipated by a shrinking, charge-leaking sphere.
    """

    # The user has provided symbolic inputs 'a' (radius) and 'V' (potential).
    # We will derive the final expression symbolically.

    print("Step 1: The total heat dissipated is equal to the initial electrostatic energy (U_initial) of the sphere.")
    print("H_total = U_initial")
    print("-" * 30)

    print("Step 2: The formula for the initial energy of a charged conductor is:")
    print("U_initial = (1/2) * Q * V")
    print("where Q is the charge and V is the potential.")
    print("-" * 30)

    print("Step 3: For a sphere of radius 'a', the charge Q can be expressed in terms of potential V:")
    print("Q = 4 * pi * epsilon_0 * a * V")
    print("where epsilon_0 is the permittivity of free space.")
    print("-" * 30)

    print("Step 4: Substitute the expression for Q into the energy formula.")
    print("U_initial = (1/2) * (4 * pi * epsilon_0 * a * V) * V")
    print("-" * 30)

    print("Step 5: Simplify the expression to get the final formula for the total heat dissipated.")
    
    # Final equation parts
    number_2_part1 = "1/2"
    number_4 = "4"
    final_number_2 = "2"
    
    print(f"H_total = ({number_2_part1}) * {number_4} * pi * epsilon_0 * a * V^2")
    print("H_total = 2 * pi * epsilon_0 * a * V^2\n")

    print("--- Final Equation Breakdown ---")
    print("The final equation for the total Joule heat (H) is: H = 2 * pi * epsilon_0 * a * V^2")
    print("Each part of the equation is:")
    print(f"Number: {final_number_2}")
    print("Constant: pi")
    print("Constant: epsilon_0 (permittivity of free space)")
    print("Variable: a (initial radius)")
    print("Variable: V (initial potential)")

# Execute the function to show the solution
solve_sphere_heat_dissipation()