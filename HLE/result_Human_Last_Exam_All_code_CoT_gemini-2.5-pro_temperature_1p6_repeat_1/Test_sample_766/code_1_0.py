import math

def find_mistake_in_crystal_data():
    """
    Calculates the density for each crystal structure dataset and compares it
    to the reported value to find the one with an error.
    """
    # Avogadro's number in mol^-1
    N_A = 6.02214076e23

    datasets = {
        'A': {'M': 608.58, 'U': 779.57, 'Z': 1, 'Dc_reported': 1.296},
        'B': {'M': 2568.09, 'U': 2618.4, 'Z': 1, 'Dc_reported': 1.629},
        'C': {'M': 1365.24, 'U': 5760.2, 'Z': 4, 'Dc_reported': 1.574},
        'D': {'M': 2804.61, 'U': 5986.6, 'Z': 2, 'Dc_reported': 1.556},
        'E': {'M': 1530.12, 'U': 5976, 'Z': 4, 'Dc_reported': 1.701}
    }

    max_error = 0
    mistake_dataset = None
    results = []

    print("Checking for inconsistencies in crystal structure data by recalculating density (Dc).\n")
    # Formula: Dc (g/cm^3) = (Z * M (g/mol)) / (U (Å^3) * 1e-24 * N_A)
    # Note: 1 g/cm^3 is numerically equal to 1 Mg/m^3

    for name, data in datasets.items():
        M = data['M']
        U = data['U']
        Z = data['Z']
        dc_reported = data['Dc_reported']
        
        # Calculate density
        dc_calculated = (Z * M) / (U * 1e-24 * N_A)
        
        # Calculate the absolute difference
        error = abs(dc_calculated - dc_reported)
        
        # Keep track of the dataset with the maximum error
        if error > max_error:
            max_error = error
            mistake_dataset = name
        
        # Store results for printing
        results.append({
            "name": name,
            "Z": Z,
            "M": M,
            "U": U,
            "dc_reported": dc_reported,
            "dc_calculated": dc_calculated
        })

    for res in results:
        print(f"Dataset {res['name']}:")
        print(f"Equation: Dc = ({res['Z']} * {res['M']}) / ({res['U']} * 1e-24 * {N_A:.6e})")
        print(f"Calculated Dc = {res['dc_calculated']:.4f} Mg m⁻³")
        print(f"Reported Dc   = {res['dc_reported']:.4f} Mg m⁻³")
        print("-" * 30)

    print(f"\nThe dataset with the largest discrepancy is '{mistake_dataset}'.")

find_mistake_in_crystal_data()
<<<E>>>