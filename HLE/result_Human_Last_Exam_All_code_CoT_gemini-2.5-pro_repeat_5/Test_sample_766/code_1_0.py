import math

def main():
    """
    Analyzes crystallographic data to find inconsistencies.
    """
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
            'beta': 105.193, 'U_reported': 5760.2
        },
        'D': {
            'system': 'monoclinic', 'a': 15.5265, 'b': 23.9138, 'c': 17.7749,
            'beta': 114.893, 'U_reported': 5986.6
        },
        'E': {
            'system': 'orthorhombic', 'a': 7.5560, 'b': 28.392, 'c': 27.854,
            'U_reported': 5976.0
        }
    }

    mistake_dataset = None
    max_diff = 0.0

    print("Verifying unit cell volume for each dataset...")
    print("-" * 60)

    for name, data in datasets.items():
        params = data
        U_calculated = 0.0

        if data['system'] == 'triclinic':
            a, b, c = params['a'], params['b'], params['c']
            alpha = math.radians(params['alpha'])
            beta = math.radians(params['beta'])
            gamma = math.radians(params['gamma'])
            cos_a, cos_b, cos_g = math.cos(alpha), math.cos(beta), math.cos(gamma)
            volume_factor = math.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g)
            U_calculated = a * b * c * volume_factor
        
        elif data['system'] == 'monoclinic':
            a, b, c = params['a'], params['b'], params['c']
            beta = math.radians(params['beta'])
            U_calculated = a * b * c * math.sin(beta)
            
        elif data['system'] == 'orthorhombic':
            a, b, c = params['a'], params['b'], params['c']
            U_calculated = a * b * c
            
        U_reported = data['U_reported']
        diff = abs(U_calculated - U_reported)
        
        print(f"Dataset {name} ({data['system']}):")
        print(f"  Reported Volume:   {U_reported:.2f} Å³")
        print(f"  Calculated Volume: {U_calculated:.2f} Å³")
        print(f"  Difference:        {diff:.2f} Å³")
        print("-" * 60)

        if diff > max_diff:
            max_diff = diff
            mistake_dataset = name

    print("\nConclusion:")
    if mistake_dataset and max_diff > 1.0:
        print(f"The largest inconsistency is found in Dataset {mistake_dataset}.")
        
        data = datasets[mistake_dataset]
        a, b, c = data['a'], data['b'], data['c']
        beta_deg = data['beta']
        beta_rad = math.radians(beta_deg)
        U_calc = a * b * c * math.sin(beta_rad)
        
        print("\nThe volume calculation for this monoclinic cell is U = a * b * c * sin(β).")
        print("Using the reported parameters:")
        print(f"U = {a} * {b} * {c} * sin({beta_deg}°)")
        print(f"U = {U_calc:.2f} Å³")
        print(f"This calculated volume ({U_calc:.2f} Å³) does not match the reported volume ({data['U_reported']:.2f} Å³).")
    else:
        print("No significant errors found in the datasets.")

if __name__ == "__main__":
    main()
    print("<<<C>>>")
