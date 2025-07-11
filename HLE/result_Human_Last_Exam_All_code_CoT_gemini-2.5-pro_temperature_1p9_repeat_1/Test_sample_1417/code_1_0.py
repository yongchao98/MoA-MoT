import pandas as pd

def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling system
    based on maximizing the radiation length to minimize particle interaction.
    """
    # Data for the materials including their radiation length (X0) in cm.
    # A longer radiation length is better as it means less particle interaction.
    materials_data = {
        'Material': ['A. Titanium', 'B. Aluminium', 'C. 316 Stainless Steel', 'D. Copper', 'E. Nickle'],
        'Radiation Length (X0) [cm]': [3.56, 8.9, 1.76, 1.43, 1.48]
    }

    df = pd.DataFrame(materials_data)
    
    # Find the material with the maximum radiation length
    optimal_material = df.loc[df['Radiation Length (X0) [cm]'].idxmax()]

    print("For a particle detector, the cooling system tubes should interact with passing particles as little as possible.")
    print("This is achieved by choosing a material with the longest possible radiation length (X0).\n")
    
    print("Comparing the radiation lengths of the given materials:")
    print(df.to_string(index=False))
    print("\n")

    print(f"The optimum choice is '{optimal_material['Material']}' because it has the highest radiation length ({optimal_material['Radiation Length (X0) [cm]']} cm),")
    print("which will minimize the unwanted scattering and energy loss of the particles being measured.")


find_optimal_material()