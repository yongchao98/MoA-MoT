import math

def calculate_volume(system, params):
    """Calculates the unit cell volume based on the crystal system and parameters."""
    if system == 'triclinic':
        a, b, c, alpha_deg, beta_deg, gamma_deg = params
        alpha = math.radians(alpha_deg)
        beta = math.radians(beta_deg)
        gamma = math.radians(gamma_deg)
        cos_a = math.cos(alpha)
        cos_b = math.cos(beta)
        cos_g = math.cos(gamma)
        
        volume_sq = 1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g
        if volume_sq < 0:
            return -1 # Should not happen for valid cell parameters
        
        volume = a * b * c * math.sqrt(volume_sq)
        return volume

    elif system == 'monoclinic':
        a, b, c, beta_deg = params
        beta = math.radians(beta_deg)
        volume = a * b * c * math.sin(beta)
        return volume

    elif system == 'orthorhombic':
        a, b, c = params
        volume = a * b * c
        return volume
    else:
        return None

# --- Data from the problem ---
datasets = {
    'A': {
        'system': 'triclinic',
        'params': (7.7810, 7.9273, 14.5543, 75.197, 88.156, 64.398),
        'given_U': 779.57
    },
    'B': {
        'system': 'triclinic',
        'params': (11.7229, 14.2639, 15.9549, 93.364, 100.301, 91.589),
        'given_U': 2618.4
    },
    'C': {
        'system': 'monoclinic',
        'params': (17.7430, 16.0855, 20.9134, 105.193),
        'given_U': 5760.2
    },
    'D': {
        'system': 'monoclinic',
        'params': (15.5265, 23.9138, 17.7749, 114.893),
        'given_U': 5986.6
    },
    'E': {
        'system': 'orthorhombic',
        'params': (7.5560, 28.392, 27.854),
        'given_U': 5976
    }
}

mistake_found_in = None

for name, data in datasets.items():
    system = data['system']
    params = data['params']
    given_U = data['given_U']
    
    calculated_U = calculate_volume(system, params)
    difference = abs(calculated_U - given_U)

    print(f"--- Checking Dataset {name} ({system}) ---")
    if system == 'triclinic':
        a, b, c, alpha, beta, gamma = params
        print(f"Formula: U = a*b*c * sqrt(1 - cos²(α) - cos²(β) - cos²(γ) + 2*cos(α)*cos(β)*cos(γ))")
        print(f"Equation: U = {a} * {b} * {c} * sqrt(1 - cos²({alpha}) - cos²({beta}) - cos²({gamma}) + 2*cos({alpha})*cos({beta})*cos({gamma}))")
    elif system == 'monoclinic':
        a, b, c, beta = params
        print(f"Formula: U = a * b * c * sin(β)")
        print(f"Equation: U = {a} * {b} * {c} * sin({beta})")
    elif system == 'orthorhombic':
        a, b, c = params
        print(f"Formula: U = a * b * c")
        print(f"Equation: U = {a} * {b} * {c}")

    print(f"Calculated U = {calculated_U:.2f} Å³")
    print(f"Given U      = {given_U:.2f} Å³")
    print(f"Difference   = {difference:.2f} Å³")
    
    # A difference greater than 0.5 is a strong indicator of an error, not just rounding.
    if difference > 0.5:
        mistake_found_in = name
        print("\n>>> A significant discrepancy is found in this dataset.")

print("\n--- Conclusion ---")
if mistake_found_in:
    print(f"The inconsistency between the calculated volume and the given volume is largest in Dataset {mistake_found_in}.")
    print("Therefore, the number that has been altered is in this dataset.")
else:
    print("No significant discrepancies were found based on the volume calculation.")
