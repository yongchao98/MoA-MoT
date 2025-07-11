import pandas as pd
from io import StringIO

def solve():
    """
    This function determines the optimal material for a particle detector's
    cooling system by finding the material with the longest radiation length,
    which minimizes interference with the particles being measured.
    """

    # Data for the materials. Radiation Length (X0) is a key parameter.
    # A longer radiation length means the material is more 'transparent' to particles.
    # Z = Atomic Number, Density is in g/cm^3, Radiation_Length is in cm.
    # For Stainless Steel, we use properties of Iron (Fe), its main component.
    data = """Material,Symbol,Z,Density (g/cm^3),Radiation_Length (cm)
Titanium,Ti,22,4.51,3.56
Aluminium,Al,13,2.70,8.9
316 Stainless Steel,Fe-based,26,8.00,1.76
Copper,Cu,29,8.96,1.43
Nickel,Ni,28,8.90,1.42
"""

    df = pd.read_csv(StringIO(data))

    print("Comparing materials based on Radiation Length (X0).")
    print("A longer radiation length is better as it means less interaction with particles.")
    print("-" * 70)
    print(df.to_string(index=False))
    print("-" * 70)

    # Find the material with the maximum radiation length
    optimal_material_row = df.loc[df['Radiation_Length (cm)'].idxmax()]
    optimal_material_name = optimal_material_row['Material']
    max_radiation_length = optimal_material_row['Radiation_Length (cm)']

    print(f"\nAnalysis:")
    print(f"The material with the longest radiation length is '{optimal_material_name}'.")
    print(f"Its radiation length is {max_radiation_length} cm, which is significantly higher than the other options.")
    print("This property makes it the most 'transparent' to passing particles, causing the least interference with detector measurements.")
    print("\nConclusion:")
    print(f"The optimum choice is {optimal_material_name}.")


solve()