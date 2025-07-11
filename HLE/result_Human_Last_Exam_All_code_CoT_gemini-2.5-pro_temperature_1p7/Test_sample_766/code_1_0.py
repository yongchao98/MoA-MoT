import math

def check_crystallography_data():
    """
    This function recalculates the density for each crystallographic dataset
    and compares it to the reported value to find inconsistencies.
    """
    
    # Avogadro's number
    N_A = 6.02214076e23
    
    # Data from the problem description
    datasets = [
        {'label': 'A', 'M': 608.58, 'Z': 1, 'U': 779.57, 'Dc_reported': 1.296},
        {'label': 'B', 'M': 2568.09, 'Z': 1, 'U': 2618.4, 'Dc_reported': 1.629},
        {'label': 'C', 'M': 1365.24, 'Z': 4, 'U': 5760.2, 'Dc_reported': 1.574},
        {'label': 'D', 'M': 2804.61, 'Z': 2, 'U': 5986.6, 'Dc_reported': 1.556},
        {'label': 'E', 'M': 1530.12, 'Z': 4, 'U': 5976, 'Dc_reported': 1.701}
    ]
    
    print("Verifying reported densities by calculation: Dc = (Z * M) / (U * N_A * 10^-24)\n")

    max_discrepancy = -1.0
    mistake_dataset = None

    for data in datasets:
        label = data['label']
        M = data['M']
        Z = data['Z']
        U = data['U']
        Dc_reported = data['Dc_reported']
        
        # Dc (g/cm^3) = (Z * M (g/mol)) / (U (Å^3) * 10^-24 (cm^3/Å^3) * N_A (mol^-1))
        # Note: 1 Mg m^-3 = 1 g cm^-3
        Dc_calculated = (Z * M) / (U * 1e-24 * N_A)
        
        discrepancy = abs(Dc_calculated - Dc_reported)
        
        # Print the full equation with the numbers used
        print(f"Dataset {label}:")
        print(f"Reported Dc = {Dc_reported} Mg m⁻³")
        print(f"Calculated Dc = ({Z} * {M}) / ({U} * {N_A:.6e} * 1e-24)")
        print(f"Result -> Dc = {Dc_calculated:.4f} Mg m⁻³")
        print(f"Discrepancy = {discrepancy:.5f}\n")
        
        if discrepancy > max_discrepancy:
            max_discrepancy = discrepancy
            mistake_dataset = label

    print("--------------------------------------------------")
    print(f"Conclusion: Dataset '{mistake_dataset}' has the largest discrepancy between its reported density and the density calculated from its other parameters.")
    print("This suggests it is the one with an altered or incorrect number.")

# Run the check
check_crystallography_data()