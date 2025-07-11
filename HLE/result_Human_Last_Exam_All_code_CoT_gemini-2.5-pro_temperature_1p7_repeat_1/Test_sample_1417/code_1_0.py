import pandas as pd

def find_optimal_material():
    """
    Analyzes materials to find the optimum choice for a particle detector cooling system.

    The unique parameter for a particle detector compared to a spacecraft is the need
    to minimize interaction with the particles being measured. This is because the
    cooling system is inside the detector and in the path of the particles.
    Interactions (like scattering or energy loss) corrupt the measurement data.

    This "transparency" to particles is best quantified by the Radiation Length (X₀).
    A longer radiation length means less interaction. Therefore, the optimum material
    is the one with the maximum radiation length.
    """
    # Data for the materials. X0 is the radiation length in cm.
    # Source: Particle Data Group (PDG)
    data = {
        'Material': ['Titanium', 'Aluminium', '316 Stainless Steel', 'Copper', 'Nickel'],
        'Choice': ['A', 'B', 'C', 'D', 'E'],
        'Primary Element Z': [22, 13, 26, 29, 28],
        'Radiation Length X0 (cm)': [3.56, 8.9, 1.76, 1.44, 1.42]
    }
    
    df = pd.DataFrame(data)
    
    # Find the material with the maximum radiation length
    optimum_material = df.loc[df['Radiation Length X0 (cm)'].idxmax()]

    print("--- Analysis of Materials for Particle Detector Cooling Tubes ---")
    print("\nUnique Parameter: Minimizing particle interaction to ensure measurement accuracy.")
    print("Physical Quantity to Maximize: Radiation Length (X0). A longer X0 means the material is more 'transparent' to particles.\n")
    
    print("Material Comparison:")
    print(df.to_string(index=False))
    
    print("\n--- Conclusion ---")
    print(f"The material with the longest radiation length is '{optimum_material['Material']}'.")
    print(f"It has a radiation length of {optimum_material['Radiation Length X0 (cm)']} cm, which is significantly higher than the other options.")
    print("This property makes it the optimum choice to minimize interference with the particle measurements.")

find_optimal_material()
<<<B>>>