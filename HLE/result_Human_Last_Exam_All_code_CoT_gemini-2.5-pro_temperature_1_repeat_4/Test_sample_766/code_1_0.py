import math

def analyze_crystal_data():
    """
    Calculates the unit cell volume from cell parameters for several datasets
    and compares it to the reported volume to find inconsistencies.
    """
    
    # Data from the problem
    datasets = {
        'A': {
            'system': 'triclinic', 'a': 7.7810, 'b': 7.9273, 'c': 14.5543,
            'alpha': 75.197, 'beta': 88.156, 'gamma': 64.398, 'U_reported': 779.57
        },
        'B': {
            'system': 'triclinic', 'a': 11.7229, 'b': 14.2639, 'c': 15.9549,
            'alpha': 93.364, 'beta': 100.301, 'gamma': 91.589, 'U_reported': 2618.4
        },
        'C': {
            'system': 'monoclinic', 'a': 17.7430, 'b': 16.0855, 'c': 20.9134,
            'beta': 105.193, 'U_reported': 5760.2, 'alpha': 90, 'gamma': 90
        },
        'D': {
            'system': 'monoclinic', 'a': 15.5265, 'b': 23.9138, 'c': 17.7749,
            'beta': 114.893, 'U_reported': 5986.6, 'alpha': 90, 'gamma': 90
        },
        'E': {
            'system': 'orthorhombic', 'a': 7.5560, 'b': 28.392, 'c': 27.854,
            'U_reported': 5976, 'alpha': 90, 'beta': 90, 'gamma': 90
        }
    }

    max_discrepancy = 0
    mistake_dataset = None

    print("Verifying the unit cell volume (U) for each dataset...")
    print("-" * 70)

    for name, data in datasets.items():
        a, b, c = data['a'], data['b'], data['c']
        alpha_deg, beta_deg, gamma_deg = data.get('alpha', 90), data.get('beta', 90), data.get('gamma', 90)
        
        alpha_rad = math.radians(alpha_deg)
        beta_rad = math.radians(beta_deg)
        gamma_rad = math.radians(gamma_deg)

        U_calc = 0
        print(f"Dataset {name} ({data['system']}):")
        
        if data['system'] == 'triclinic':
            # Formula: U = a*b*c * sqrt(1 - cos²α - cos²β - cos²γ + 2*cosα*cosβ*cosγ)
            cos_a = math.cos(alpha_rad)
            cos_b = math.cos(beta_rad)
            cos_g = math.cos(gamma_rad)
            sqrt_term = math.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g)
            U_calc = a * b * c * sqrt_term
            print(f"U = {a} * {b} * {c} * sqrt(1 - cos²({alpha_deg}) - cos²({beta_deg}) - cos²({gamma_deg}) + 2*cos({alpha_deg})*cos({beta_deg})*cos({gamma_deg}))")

        elif data['system'] == 'monoclinic':
            # Formula: U = a*b*c * sinβ
            U_calc = a * b * c * math.sin(beta_rad)
            print(f"U = {a} * {b} * {c} * sin({beta_deg})")
            
        elif data['system'] == 'orthorhombic':
            # Formula: U = a*b*c
            U_calc = a * b * c
            print(f"U = {a} * {b} * {c}")

        U_reported = data['U_reported']
        discrepancy = abs(U_calc - U_reported)
        
        print(f"Calculated U = {U_calc:.2f} Å³")
        print(f"Reported U   = {U_reported:.2f} Å³")
        print(f"Discrepancy  = {discrepancy:.2f} Å³\n")
        
        if discrepancy > max_discrepancy:
            max_discrepancy = discrepancy
            mistake_dataset = name

    print("-" * 70)
    print(f"The dataset with the largest discrepancy is '{mistake_dataset}'.")
    print(f"The difference between its reported and calculated volume is {max_discrepancy:.2f} Å³.")

analyze_crystal_data()