import math

def find_mistake():
    """
    Calculates the density for each crystal structure dataset and compares it
    to the given value to find the dataset with an error.
    """
    # Avogadro's constant (mol⁻¹)
    N_A = 6.02214076e23

    datasets = {
        'A': {'M': 608.58, 'U': 779.57, 'Z': 1, 'Dc_given': 1.296},
        'B': {'M': 2568.09, 'U': 2618.4, 'Z': 1, 'Dc_given': 1.629},
        'C': {'M': 1365.24, 'U': 5760.2, 'Z': 4, 'Dc_given': 1.574},
        'D': {'M': 2804.61, 'U': 5986.6, 'Z': 2, 'Dc_given': 1.556},
        'E': {'M': 1530.12, 'U': 5976, 'Z': 4, 'Dc_given': 1.701}
    }

    max_relative_error = 0
    error_dataset_name = None

    print("--- Verifying Calculated Density (Dc) for Each Dataset ---")

    for name, data in datasets.items():
        M = data['M']
        U = data['U']
        Z = data['Z']
        Dc_given = data['Dc_given']

        # Calculate Dc in g/cm³ (equivalent to Mg/m³)
        # Dc = (Z * M) / (U * N_A * 1e-24)
        Dc_calc = (Z * M) / (U * N_A * 1e-24)
        
        relative_error = abs((Dc_calc - Dc_given) / Dc_given)

        print(f"\nDataset {name}:")
        print(f"  Given Dc    = {Dc_given:.4f} Mg m⁻³")
        print(f"  Calculated Dc = {Dc_calc:.4f} Mg m⁻³")
        print(f"  Relative Error = {relative_error:.4%}")

        if relative_error > max_relative_error:
            max_relative_error = relative_error
            error_dataset_name = name

    print("\n--- Conclusion ---")
    print(f"The dataset with the largest discrepancy is '{error_dataset_name}'.")
    print("This indicates a likely error or altered number in this dataset's parameters.")
    
    # Show the calculation for the identified dataset
    error_data = datasets[error_dataset_name]
    M = error_data['M']
    U = error_data['U']
    Z = error_data['Z']
    Dc_calc = (Z * M) / (U * N_A * 1e-24)

    print("\nHere is the calculation for the incorrect dataset:")
    print(f"Dc = (Z * M) / (U * N_A * 10⁻²⁴)")
    print(f"Calculated Dc = ({Z} * {M}) / ({U} * {N_A} * 1e-24) = {Dc_calc:.4f}")
    print(f"The given Dc is {error_data['Dc_given']:.4f}, which does not match the calculated value as closely as in other datasets.")
    
find_mistake()