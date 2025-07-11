import pandas as pd

def find_optimal_material():
    """
    This function compares different materials based on their radiation length
    to find the optimal choice for a particle detector cooling system.
    """
    # Data for the materials. Radiation Length (X0) is a key parameter.
    # A longer radiation length is better as it means less interaction with particles.
    # Data sourced from the Particle Data Group (PDG).
    # For 316 Stainless Steel, properties are approximated by its main component, Iron (Fe).
    data = {
        'Material': ['Titanium', 'Aluminium', '316 Stainless Steel (Fe)', 'Copper', 'Nickel'],
        'Symbol': ['Ti', 'Al', 'Fe', 'Cu', 'Ni'],
        'Atomic Number (Z)': [22, 13, 26, 29, 28],
        'Density (g/cm^3)': [4.51, 2.70, 7.87, 8.96, 8.91],
        'Radiation Length X0 (cm)': [3.60, 8.90, 1.76, 1.43, 1.42]
    }

    df = pd.DataFrame(data)

    print("Comparing materials for particle detector cooling tubes...")
    print("The unique and most important parameter is maximizing the radiation length (X0)")
    print("to minimize interaction with the particles being measured.\n")
    print(df.to_string(index=False))

    # Find the material with the maximum radiation length
    optimal_material_row = df.loc[df['Radiation Length X0 (cm)'].idxmax()]
    optimal_material_name = optimal_material_row['Material'].split(' ')[0] # Get 'Aluminium' from 'Aluminium'

    print(f"\n\nConclusion:")
    print(f"Comparing the radiation lengths, {optimal_material_name} has the highest value ({optimal_material_row['Radiation Length X0 (cm)']} cm).")
    print(f"This means it will interfere the least with the particle measurements, making it the optimum choice.")

find_optimal_material()