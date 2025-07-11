import math

def check_crystal_data():
    """
    Checks the consistency of crystallographic data by recalculating the unit cell volume
    from the lattice parameters and comparing it to the given volume.
    """

    # Data from the problem description
    datasets = {
        'A': {
            "system": "triclinic", "params": (7.7810, 7.9273, 14.5543, 75.197, 88.156, 64.398), 
            "U_given": 779.57
        },
        'B': {
            "system": "triclinic", "params": (11.7229, 14.2639, 15.9549, 93.364, 100.301, 91.589), 
            "U_given": 2618.4
        },
        'C': {
            "system": "monoclinic", "params": (17.7430, 16.0855, 20.9134, 90, 105.193, 90),
            "U_given": 5760.2
        },
        'D': {
            "system": "monoclinic", "params": (15.5265, 23.9138, 17.7749, 90, 114.893, 90),
            "U_given": 5986.6
        },
        'E': {
            "system": "orthorhombic", "params": (7.5560, 28.392, 27.854, 90, 90, 90),
            "U_given": 5976
        },
    }

    mistake_found_in = None
    max_diff = -1

    for name, data in datasets.items():
        system = data["system"]
        params = data["params"]
        a, b, c, alpha_deg, beta_deg, gamma_deg = params
        U_given = data["U_given"]
        
        U_calc = 0
        
        print(f"--- Analyzing Dataset {name} ({system}) ---")
        print(f"Given Volume U = {U_given} Å³")
        
        if system == "triclinic":
            alpha_rad = math.radians(alpha_deg)
            beta_rad = math.radians(beta_deg)
            gamma_rad = math.radians(gamma_deg)
            
            cos_a = math.cos(alpha_rad)
            cos_b = math.cos(beta_rad)
            cos_g = math.cos(gamma_rad)

            sqrt_term = math.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g)
            U_calc = a * b * c * sqrt_term
            
            print("Calculation: U = a * b * c * sqrt(1 - cos²(α) - cos²(β) - cos²(γ) + 2*cos(α)*cos(β)*cos(γ))")
            print(f"U_calc = {a} * {b} * {c} * sqrt(1 - cos²({alpha_deg}) - cos²({beta_deg}) - cos²({gamma_deg}) + 2*cos({alpha_deg})*cos({beta_deg})*cos({gamma_deg}))")
        
        elif system == "monoclinic":
            beta_rad = math.radians(beta_deg)
            U_calc = a * b * c * math.sin(beta_rad)
            print("Calculation: U = a * b * c * sin(β)")
            print(f"U_calc = {a} * {b} * {c} * sin({beta_deg})")
        
        elif system == "orthorhombic":
            U_calc = a * b * c
            print("Calculation: U = a * b * c")
            print(f"U_calc = {a} * {b} * {c}")

        difference = abs(U_given - U_calc)
        print(f"Calculated Volume U = {U_calc:.2f} Å³")
        print(f"Difference = |{U_given} - {U_calc:.2f}| = {difference:.2f}")

        if difference > max_diff:
            max_diff = difference
            # Setting a threshold to consider it a mistake, e.g., > 0.1
            if max_diff > 0.1:
                mistake_found_in = name
        
        print("-" * (len(name) + 29))
        print()

    print(f"\nConclusion: The largest discrepancy between the given and calculated volume is in Dataset {mistake_found_in}.")
    print("This indicates an error in the reported numbers for this dataset.")

if __name__ == '__main__':
    check_crystal_data()
<<<A>>>