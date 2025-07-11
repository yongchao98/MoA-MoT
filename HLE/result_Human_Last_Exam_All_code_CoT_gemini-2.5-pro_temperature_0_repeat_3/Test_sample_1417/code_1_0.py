import pandas as pd

def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling system
    by maximizing the radiation length.
    """
    # Data for the materials. Radiation Length (X0) is the key parameter.
    # Density is in g/cm^3, and Radiation Length is in cm.
    # Stainless steel is an alloy, but its properties are dominated by Iron.
    data = {
        'Material': ['Titanium', 'Aluminium', '316 Stainless Steel', 'Copper', 'Nickel'],
        'Symbol': ['Ti', 'Al', 'Fe-based alloy', 'Cu', 'Ni'],
        'Atomic Number (Z)': [22, 13, '~26 (Iron)', 29, 28],
        'Density (g/cm^3)': [4.51, 2.70, 8.0, 8.96, 8.91],
        'Radiation Length (cm)': [3.6, 8.9, 1.76, 1.43, 1.42]
    }

    df = pd.DataFrame(data)

    print("Comparing materials based on Radiation Length (X0):")
    print("A longer radiation length is better as it means less interaction with particles.")
    print("-" * 60)
    print(df.to_string(index=False))
    print("-" * 60)

    # Find the material with the maximum radiation length
    optimal_material_row = df.loc[df['Radiation Length (cm)'].idxmax()]
    optimal_material_name = optimal_material_row['Material']
    optimal_value = optimal_material_row['Radiation Length (cm)']

    print(f"\nThe material with the maximum radiation length is '{optimal_material_name}' with a value of {optimal_value} cm.")
    print("Therefore, Aluminium is the optimum choice to minimize interference with particle measurements.")

find_optimal_material()