import math

def find_optimal_material():
    """
    Calculates the radiation length in cm for various materials to determine the
    optimal choice for a particle detector cooling system, where minimizing particle
    interaction is the key parameter.
    """
    # Material properties:
    # Name, Density (g/cm^3), Radiation Length (g/cm^2)
    # Data sourced from the Particle Data Group (PDG)
    materials = [
        {'name': 'Titanium', 'density': 4.51, 'rad_length_mass': 16.0},
        {'name': 'Aluminium', 'density': 2.70, 'rad_length_mass': 24.01},
        {'name': '316 Stainless Steel', 'density': 8.00, 'rad_length_mass': 13.9}, # Approximated as Iron/Steel
        {'name': 'Copper', 'density': 8.96, 'rad_length_mass': 12.86},
        {'name': 'Nickel', 'density': 8.91, 'rad_length_mass': 12.96}
    ]

    print("The unique parameter to maximize for a particle detector cooling system is the radiation length (X0).")
    print("A longer radiation length minimizes the interaction of particles with the material.")
    print("Calculating the radiation length in cm for pipes of similar thickness...\n")
    print("Formula: Radiation Length [cm] = Radiation Length [g/cm^2] / Density [g/cm^3]\n")

    best_material = None
    max_rad_length_cm = -1.0

    for material in materials:
        name = material['name']
        density = material['density']
        rad_length_mass = material['rad_length_mass']
        
        # Calculate radiation length in cm
        rad_length_cm = rad_length_mass / density
        
        print(f"For {name}:")
        # Print the equation with all numbers
        print(f"  Radiation Length [cm] = {rad_length_mass} / {density} = {rad_length_cm:.2f} cm")
        
        if rad_length_cm > max_rad_length_cm:
            max_rad_length_cm = rad_length_cm
            best_material = name

    print(f"\nConclusion: {best_material} has the longest radiation length ({max_rad_length_cm:.2f} cm) and is therefore the optimum choice.")

if __name__ == '__main__':
    find_optimal_material()