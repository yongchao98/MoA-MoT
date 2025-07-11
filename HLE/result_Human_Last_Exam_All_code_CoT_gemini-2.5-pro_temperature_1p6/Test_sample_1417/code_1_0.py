import sys

def solve_material_choice():
    """
    This script determines the optimum material for a particle detector's
    integrated cooling system by focusing on the unique parameter relevant
    to this application.
    """

    # Step 1: Explain the unique parameter.
    print("For a particle detector, the cooling system is part of the 'material budget'.")
    print("It is critical to choose a material that minimizes interactions with the particles being detected.")
    print("This property is quantified by the Radiation Length (X0). A longer radiation length is better.")
    print("The goal is to find the material with the maximum radiation length.\n")

    # Step 2: Define material properties.
    # Data includes Atomic Number (Z) and Radiation Length (X0) in cm.
    # 316 Stainless Steel is an alloy, mostly Iron (Fe), so we use Iron's properties as a close approximation.
    materials = [
        {'name': 'Aluminium', 'symbol': 'Al', 'Z': 13, 'X0_cm': 8.897},
        {'name': 'Titanium', 'symbol': 'Ti', 'Z': 22, 'X0_cm': 3.560},
        {'name': '316 Stainless Steel (approx. Iron)', 'symbol': 'Fe', 'Z': 26, 'X0_cm': 1.757},
        {'name': 'Copper', 'symbol': 'Cu', 'Z': 29, 'X0_cm': 1.436},
        {'name': 'Nickel', 'symbol': 'Ni', 'Z': 28, 'X0_cm': 1.424}
    ]

    # Step 3: Find the material with the maximum radiation length.
    optimal_material = None
    max_x0 = -1

    print("--- Material Properties ---")
    for material in materials:
        print(f"Material: {material['name']:<35} | Atomic Number (Z): {material['Z']:<5} | Radiation Length (X0): {material['X0_cm']:.3f} cm")
        if material['X0_cm'] > max_x0:
            max_x0 = material['X0_cm']
            optimal_material = material
    print("---------------------------\n")


    # Step 4: Print the conclusion.
    if optimal_material:
        print("Conclusion:")
        print(f"Comparing the radiation lengths, {optimal_material['name']} has the highest value ({optimal_material['X0_cm']:.3f} cm).")
        print("This means it will interfere the least with the particles being measured.")
        print(f"Therefore, {optimal_material['name']} is the optimum choice for this unique requirement.")
    else:
        print("Could not determine the optimal material.", file=sys.stderr)

solve_material_choice()

# Final Answer choice B corresponds to Aluminium.
<<<B>>>