import operator

def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling system
    by finding the material with the maximum radiation length.
    """
    # Data: Material and its approximate Radiation Length (X₀) in cm.
    # A longer radiation length is better as it implies less particle interaction.
    materials = {
        "Titanium": 3.56,
        "Aluminium": 8.9,
        "316 Stainless Steel": 1.76, # Approximate value for steel
        "Copper": 1.43,
        "Nickel": 1.48
    }

    print("The goal is to maximize the radiation length (X₀) to minimize particle interaction.")
    print("Comparing the radiation lengths of the candidate materials:\n")

    # Print the values for each material
    for material, x0 in materials.items():
        print(f"{material:<20}: X₀ = {x0} cm")

    # Find the material with the maximum radiation length
    optimal_material = max(materials.items(), key=operator.itemgetter(1))

    print(f"\nBased on the data, {optimal_material[0]} has the highest radiation length.")
    print(f"Therefore, it is the optimum choice to minimize the impact on particle measurements.")

find_optimal_material()