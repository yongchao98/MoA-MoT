import pandas as pd
from io import StringIO

def find_optimal_material():
    """
    Compares materials for a particle detector's cooling system based on the unique
    parameter of minimizing particle interaction. The best material will have a low
    atomic number (Z) and a long radiation length (X0).
    """

    # Data for the materials
    # Z = Atomic Number
    # Density is in g/cm^3
    # X0 (Radiation Length) is in cm
    data = """Material,Z_approx,Density_g_cm3,Radiation_Length_X0_cm
Aluminium,13,2.70,8.9
Titanium,22,4.51,3.6
316 Stainless Steel,26,8.00,1.76
Nickel,28,8.91,1.42
Copper,29,8.96,1.43
"""

    # Create a pandas DataFrame
    df = pd.read_csv(StringIO(data))

    # Find the material with the maximum radiation length
    optimal_material = df.loc[df['Radiation_Length_X0_cm'].idxmax()]

    # Print the comparison table
    print("--- Material Properties for Particle Detector Cooling System ---")
    print(df.to_string(index=False))
    print("\n--- Analysis ---")
    print("The unique parameter for a particle detector is minimizing interaction with the particles being measured.")
    print("This requires a material with a low Atomic Number (Z) and a long Radiation Length (X0).")
    print("\nBased on the data, the optimal choice is the material with the longest radiation length.")
    print(f"\nOptimal Choice: {optimal_material['Material']}")
    print(f"It has the longest Radiation Length of {optimal_material['Radiation_Length_X0_cm']} cm and the lowest Atomic Number.")

if __name__ == '__main__':
    find_optimal_material()