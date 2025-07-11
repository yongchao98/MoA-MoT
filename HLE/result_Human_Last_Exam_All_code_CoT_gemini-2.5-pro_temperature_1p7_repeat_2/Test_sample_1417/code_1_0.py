import operator

def find_optimal_material():
    """
    Calculates the radiation length for various materials to determine the best choice
    for a particle detector cooling system, where minimizing particle interaction is key.
    Data is sourced from the Particle Data Group (PDG).
    """

    # Material properties:
    # Density (rho) in g/cm^3
    # Radiation Mass (X0 * rho) in g/cm^2
    materials = {
        'Titanium':          {'rho': 4.54, 'rad_mass': 16.2},
        'Aluminium':         {'rho': 2.70, 'rad_mass': 24.01},
        '316 Stainless Steel': {'rho': 8.00, 'rad_mass': 13.9}, # Using properties for Iron as a close approximation
        'Copper':            {'rho': 8.96, 'rad_mass': 12.86},
        'Nickle':            {'rho': 8.90, 'rad_mass': 12.6}
    }

    print("Calculating Radiation Length (X0) for each material...")
    print("The optimal material is the one with the highest radiation length.\n")

    results = {}
    for name, properties in materials.items():
        rho = properties['rho']
        rad_mass = properties['rad_mass']
        
        # Radiation Length (cm) = Radiation Mass (g/cm^2) / Density (g/cm^3)
        rad_length = rad_mass / rho
        results[name] = rad_length
        
        # The final equation is: Radiation Length = Radiation Mass / Density
        print(f"Material: {name}")
        print(f"Equation: {rad_length:.2f} cm = {rad_mass} g/cm^2 / {rho} g/cm^3")
        print("-" * 30)

    # Find the material with the maximum radiation length
    optimal_material = max(results.items(), key=operator.itemgetter(1))

    print(f"\nConclusion:")
    print(f"The material with the longest radiation length is '{optimal_material[0]}' ({optimal_material[1]:.2f} cm).")
    print("Therefore, Aluminium is the optimum choice to minimize interaction with particles being detected.")

if __name__ == "__main__":
    find_optimal_material()
<<<B>>>