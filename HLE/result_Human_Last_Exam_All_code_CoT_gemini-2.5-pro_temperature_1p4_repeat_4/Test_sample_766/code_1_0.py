import math

def check_crystallographic_data():
    """
    Calculates the density for each crystal structure dataset and compares it to the reported value.
    The dataset with the largest discrepancy is identified as having an error.
    """

    # Avogadro's constant * 10^-24 (for converting Å³ to cm³ and g/mol to g)
    NA_CONST = 6.02214076e23 * 1e-24

    datasets = [
        {'id': 'A', 'M': 608.58, 'Z': 1, 'U': 779.57, 'Dc_rep': 1.296},
        {'id': 'B', 'M': 2568.09, 'Z': 1, 'U': 2618.4, 'Dc_rep': 1.629},
        {'id': 'C', 'M': 1365.24, 'Z': 4, 'U': 5760.2, 'Dc_rep': 1.574},
        {'id': 'D', 'M': 2804.61, 'Z': 2, 'U': 5986.6, 'Dc_rep': 1.556},
        {'id': 'E', 'M': 1530.12, 'Z': 4, 'U': 5976, 'Dc_rep': 1.701},
    ]

    print("Checking consistency of crystallographic data by recalculating density (Dc).")
    print("Formula: Dc = (Z * M) / (U * N_A_const)")
    print(f"Using N_A_const = {NA_CONST:.6f}\n")
    
    max_diff = 0
    mistake_id = None

    for data in datasets:
        # Unpack data for clarity
        d_id = data['id']
        Z = data['Z']
        M = data['M']
        U = data['U']
        Dc_rep = data['Dc_rep']

        # Calculate density
        Dc_calc = (Z * M) / (U * NA_CONST)
        
        # Calculate the absolute difference
        diff = abs(Dc_calc - Dc_rep)
        
        # Keep track of the one with the largest difference
        if diff > max_diff:
            max_diff = diff
            mistake_id = d_id
            
        print(f"Dataset {d_id}:")
        print(f"Equation: ({Z} * {M}) / ({U} * {NA_CONST:.6f}) = {Dc_calc:.4f}")
        print(f"Calculated Dc = {Dc_calc:.4f} Mg m⁻³")
        print(f"Reported Dc   = {Dc_rep:.4f} Mg m⁻³")
        print(f"Difference      = {diff:.4f}\n")

    print("="*30)
    print(f"Conclusion:")
    print(f"The largest discrepancy is found in Dataset {mistake_id}.")
    print("This indicates that one of the reported numbers in this dataset is incorrect.")
    print("="*30)

check_crystallographic_data()
<<<E>>>