import operator

def solve():
    """
    This function determines the optimal material for a particle detector's
    cooling system by maximizing the radiation length (X₀), which minimizes
    particle interaction.
    """

    # Data: Radiation length (X₀) in cm for each material.
    # A longer radiation length is better as it signifies less interaction.
    # X₀ (cm) = X₀ (g/cm²) / density (g/cm³)
    materials = {
        'Aluminium': 8.9,          # Z=13, relatively low density
        'Titanium': 3.6,           # Z=22
        '316 Stainless Steel': 1.76, # Mostly Iron (Z=26)
        'Copper': 1.43,            # Z=29
        'Nickle': 1.42             # Z=28
    }

    # Find the material with the maximum radiation length
    best_material, max_x0 = max(materials.items(), key=operator.itemgetter(1))

    # --- Output ---
    print("The unique requirement for a particle detector is to minimize interaction with passing particles.")
    print("This is achieved by selecting a material with the longest possible radiation length (X₀).\n")
    print("Comparing the radiation lengths of the choices:")

    for material, x0 in materials.items():
        print(f"- {material}: X₀ = {x0:.2f} cm")

    print(f"\nBased on this unique parameter, {best_material} is the optimum choice because it has the highest radiation length.")
    print(f"\nFinal Equation: {best_material} X₀({max_x0}) = max(Aluminium X₀({materials['Aluminium']}), Titanium X₀({materials['Titanium']}), 316 Stainless Steel X₀({materials['316 Stainless Steel']}), Copper X₀({materials['Copper']}), Nickle X₀({materials['Nickle']}))")


solve()
<<<B>>>