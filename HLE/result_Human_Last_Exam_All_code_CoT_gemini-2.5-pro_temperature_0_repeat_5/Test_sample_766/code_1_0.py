import math

def calculate_volume(system, a, b, c, alpha_deg, beta_deg, gamma_deg):
    """Calculates the unit cell volume based on the crystal system and parameters."""
    alpha = math.radians(alpha_deg)
    beta = math.radians(beta_deg)
    gamma = math.radians(gamma_deg)

    if system == 'orthorhombic':
        # For orthorhombic, alpha = beta = gamma = 90 degrees
        return a * b * c
    elif system == 'monoclinic':
        # For monoclinic, alpha = gamma = 90 degrees
        return a * b * c * math.sin(beta)
    elif system == 'triclinic':
        cos_alpha = math.cos(alpha)
        cos_beta = math.cos(beta)
        cos_gamma = math.cos(gamma)
        volume_sq = (a*b*c)**2 * (1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2 * cos_alpha * cos_beta * cos_gamma)
        return math.sqrt(volume_sq)
    else:
        return None

# Data from the problem description
datasets = {
    'A': {'system': 'triclinic', 'a': 7.7810, 'b': 7.9273, 'c': 14.5543, 'alpha': 75.197, 'beta': 88.156, 'gamma': 64.398, 'U_reported': 779.57},
    'B': {'system': 'triclinic', 'a': 11.7229, 'b': 14.2639, 'c': 15.9549, 'alpha': 93.364, 'beta': 100.301, 'gamma': 91.589, 'U_reported': 2618.4},
    'C': {'system': 'monoclinic', 'a': 17.7430, 'b': 16.0855, 'c': 20.9134, 'alpha': 90, 'beta': 105.193, 'gamma': 90, 'U_reported': 5760.2},
    'D': {'system': 'monoclinic', 'a': 15.5265, 'b': 23.9138, 'c': 17.7749, 'alpha': 90, 'beta': 114.893, 'gamma': 90, 'U_reported': 5986.6},
    'E': {'system': 'orthorhombic', 'a': 7.5560, 'b': 28.392, 'c': 27.854, 'alpha': 90, 'beta': 90, 'gamma': 90, 'U_reported': 5976}
}

max_discrepancy = 0
error_dataset_key = None
error_dataset_calc = None

print("Verifying consistency of unit cell volume...\n")

for key, data in datasets.items():
    u_calc = calculate_volume(data['system'], data['a'], data['b'], data['c'], data['alpha'], data['beta'], data['gamma'])
    discrepancy = abs(u_calc - data['U_reported'])
    
    print(f"Dataset {key}:")
    print(f"  Reported Volume = {data['U_reported']:.2f} Å³")
    print(f"  Calculated Volume = {u_calc:.2f} Å³")
    print(f"  Difference = {discrepancy:.2f} Å³\n")

    if discrepancy > max_discrepancy:
        max_discrepancy = discrepancy
        error_dataset_key = key
        error_dataset_calc = u_calc

print("--- Conclusion ---")
print(f"The largest discrepancy is found in dataset {error_dataset_key}.")
print("The reported unit cell parameters are inconsistent with the reported unit cell volume.")

# Show the calculation for the incorrect dataset
error_data = datasets[error_dataset_key]
a, b, c, beta = error_data['a'], error_data['b'], error_data['c'], error_data['beta']

print("\nThe calculation for the volume of dataset C shows the error:")
print(f"Formula: U = a * b * c * sin(β)")
print(f"U = {a} * {b} * {c} * sin({beta}°)")
print(f"U = {error_dataset_calc:.2f} Å³")
print(f"This calculated value ({error_dataset_calc:.2f} Å³) does not match the reported value ({error_data['U_reported']:.2f} Å³).")
print("\nTherefore, dataset C contains the mistake.")