import math

def calculate_density(M, Z, U_A3):
    """
    Calculates the density in Mg/m^3 (which is numerically equal to g/cm^3).
    Dc = (Z * M) / (N_A * U)
    - Z: number of formula units in the unit cell
    - M: molar mass (g/mol)
    - N_A: Avogadro's number (mol^-1)
    - U: unit cell volume (Å^3)
    We need to convert U from Å^3 to cm^3 by multiplying by 1e-24.
    """
    N_A = 6.02214076e23
    density = (Z * M) / (N_A * U_A3 * 1e-24)
    return density

# Data from the problem
datasets = {
    'A': {'M': 608.58, 'Z': 1, 'U': 779.57, 'Dc_reported': 1.296},
    'B': {'M': 2568.09, 'Z': 1, 'U': 2618.4, 'Dc_reported': 1.629},
    'C': {'M': 1365.24, 'Z': 4, 'U': 5760.2, 'Dc_reported': 1.574},
    'D': {'M': 2804.61, 'Z': 2, 'U': 5986.6, 'Dc_reported': 1.556},
    'E': {'M': 1530.12, 'Z': 4, 'U': 5976, 'Dc_reported': 1.701}
}

mistake_found_in = None

# Perform calculations for each dataset
for name, data in datasets.items():
    M = data['M']
    Z = data['Z']
    U = data['U']
    Dc_reported = data['Dc_reported']
    
    # Calculate density
    Dc_calculated = calculate_density(M, Z, U)
    
    print(f"--- Dataset {name} ---")
    print(f"Reported Density: {Dc_reported} Mg m⁻³")
    # The final equation with numbers is printed here
    print(f"Calculated Density = ({Z} * {M}) / (N_A * {U} Å³)")
    print(f"                   = {Dc_calculated:.4f} Mg m⁻³")
    
    # Check for discrepancy (allowing for small rounding differences)
    if not math.isclose(Dc_calculated, Dc_reported, rel_tol=0.001): # 0.1% tolerance
        print(f"*** Discrepancy found in Dataset {name}! ***")
        mistake_found_in = name
    
    print("-" * 20 + "\n")

if mistake_found_in:
    print(f"Conclusion: The mistake is in Dataset {mistake_found_in}.")
    print("The calculated density does not match the reported value, indicating an error in one of the reported numbers for this dataset.")
else:
    print("Conclusion: No significant discrepancies found based on density calculation.")
