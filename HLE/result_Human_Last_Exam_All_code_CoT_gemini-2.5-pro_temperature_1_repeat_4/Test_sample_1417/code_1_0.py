import operator

def solve():
    """
    Analyzes materials for a particle detector's cooling system to find the
    optimum choice based on the unique parameter of radiation length.
    """
    
    # Data for materials: Name, Density (g/cm^3), Radiation Length (cm)
    # The unique parameter for a particle detector is the need to minimize interaction
    # with the particles being measured. This is achieved by maximizing the
    # radiation length (X0). A longer radiation length means the material is more
    # "transparent" to particles.
    materials = {
        'Titanium': {'density': 4.54, 'radiation_length': 3.56},
        'Aluminium': {'density': 2.70, 'radiation_length': 8.9},
        '316 Stainless Steel': {'density': 8.00, 'radiation_length': 1.74}, # Approx. for iron
        'Copper': {'density': 8.96, 'radiation_length': 1.43},
        'Nickel': {'density': 8.90, 'radiation_length': 1.45}
    }

    print("Comparing materials based on Radiation Length (X0).")
    print("A longer radiation length is better as it minimizes interference with particle measurements.\n")
    
    # Find the material with the maximum radiation length
    best_material = max(materials.items(), key=lambda item: item[1]['radiation_length'])

    for name, properties in materials.items():
        print(f"Material: {name}")
        print(f"  - Density: {properties['density']:.2f} g/cm^3")
        print(f"  - Radiation Length: {properties['radiation_length']:.2f} cm")
        if name == best_material[0]:
            print("  - <<< This material has the longest radiation length.")
        print("-" * 20)

    print(f"\nConclusion:")
    print(f"The unique parameter for a particle detector cooling system is the 'material budget',")
    print(f"which must be minimized to avoid disrupting particle measurements.")
    print(f"This is achieved by selecting the material with the maximum radiation length.")
    print(f"Based on the data, {best_material[0]} has the longest radiation length ({best_material[1]['radiation_length']:.2f} cm) and is the optimum choice.")

solve()
<<<B>>>