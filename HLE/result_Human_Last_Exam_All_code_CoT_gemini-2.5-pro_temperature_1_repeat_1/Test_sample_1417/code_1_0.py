import math

def solve_material_choice():
    """
    Calculates the radiation length for different materials to determine the best choice
    for a particle detector cooling system.
    The goal is to maximize the radiation length to minimize particle interaction.
    """
    # Material properties: Name, X0 (g/cm^2), Density (g/cm^3)
    # X0 (radiation length in g/cm^2) values are from the Particle Data Group (PDG)
    materials = {
        'Titanium': {'x0_density': 16.0, 'density': 4.51},
        'Aluminium': {'x0_density': 24.01, 'density': 2.70},
        '316 Stainless Steel': {'x0_density': 13.84, 'density': 8.00}, # Using Iron as a proxy for X0, typical density for 316 SS
        'Copper': {'x0_density': 12.86, 'density': 8.96},
        'Nickel': {'x0_density': 12.6, 'density': 8.91}
    }

    print("The unique parameter for a particle detector is minimizing particle interaction, which means maximizing the radiation length (L_rad).")
    print("The formula is: L_rad (cm) = X0 (g/cm^2) / Density (g/cm^3)\n")

    best_material = None
    max_l_rad = -1

    for name, props in materials.items():
        x0_density = props['x0_density']
        density = props['density']
        l_rad = x0_density / density
        
        # Output the calculation for each material
        print(f"Calculating for {name}:")
        print(f"L_rad = {x0_density} / {density} = {l_rad:.2f} cm")
        print("-" * 20)

        if l_rad > max_l_rad:
            max_l_rad = l_rad
            best_material = name

    print(f"\nThe material with the highest radiation length is {best_material} with L_rad = {max_l_rad:.2f} cm.")
    print("Therefore, Aluminium is the optimum choice.")

solve_material_choice()