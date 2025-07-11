import pandas as pd

def find_optimal_material():
    """
    Analyzes materials to find the one with the longest radiation length,
    which is the optimal choice for minimizing particle interaction in a detector.
    """
    # Data for the materials. Radiation Length (X0) is the key parameter.
    # A longer radiation length is better as it signifies less interaction.
    # Data is for common alloys and elements.
    # Source: Particle Data Group (pdg.lbl.gov)
    materials_data = {
        'A. Titanium': {'Radiation Length (cm)': 3.56},
        'B. Aluminium': {'Radiation Length (cm)': 8.9},
        'C. 316 Stainless Steel': {'Radiation Length (cm)': 1.76}, # Approximate value for steel
        'D. Copper': {'Radiation Length (cm)': 1.43},
        'E. Nickle': {'Radiation Length (cm)': 1.42}
    }

    print("Comparing materials based on Radiation Length (X0).")
    print("A longer radiation length is desirable to minimize particle interaction.\n")

    # Using pandas for a clean output table
    df = pd.DataFrame.from_dict(materials_data, orient='index')
    print(df)
    
    # Find the material with the maximum radiation length
    optimal_material_name = None
    max_radiation_length = -1

    for name, properties in materials_data.items():
        if properties['Radiation Length (cm)'] > max_radiation_length:
            max_radiation_length = properties['Radiation Length (cm)']
            optimal_material_name = name

    print(f"\n--- Analysis ---")
    print(f"The unique requirement for a particle detector is to minimize interaction with the particles being measured.")
    print(f"This is achieved by selecting the material with the maximum radiation length.")
    print(f"The material with the maximum radiation length is '{optimal_material_name}' with a value of {max_radiation_length} cm.")
    print("\nTherefore, Aluminium is the optimum choice.")

find_optimal_material()