import math

def check_crystal_data():
    """
    Recalculates the unit cell volume for each dataset and compares it to the reported value.
    Identifies the dataset with a significant discrepancy.
    """
    
    datasets = [
        {
            "name": "A",
            "system": "triclinic",
            "a": 7.7810, "b": 7.9273, "c": 14.5543,
            "alpha": 75.197, "beta": 88.156, "gamma": 64.398,
            "U_reported": 779.57
        },
        {
            "name": "B",
            "system": "triclinic",
            "a": 11.7229, "b": 14.2639, "c": 15.9549,
            "alpha": 93.364, "beta": 100.301, "gamma": 91.589,
            "U_reported": 2618.4
        },
        {
            "name": "C",
            "system": "monoclinic",
            "a": 17.7430, "b": 16.0855, "c": 20.9134,
            "beta": 105.193,
            "U_reported": 5760.2
        },
        {
            "name": "D",
            "system": "monoclinic",
            "a": 15.5265, "b": 23.9138, "c": 17.7749,
            "beta": 114.893,
            "U_reported": 5986.6
        },
        {
            "name": "E",
            "system": "orthorhombic",
            "a": 7.5560, "b": 28.392, "c": 27.854,
            "U_reported": 5976.0
        }
    ]

    for data in datasets:
        name = data["name"]
        system = data["system"]
        U_reported = data["U_reported"]
        
        print(f"--- Checking Dataset {name} ({system}) ---")

        if system == "triclinic":
            a, b, c = data["a"], data["b"], data["c"]
            alpha_deg, beta_deg, gamma_deg = data["alpha"], data["beta"], data["gamma"]
            
            alpha_rad = math.radians(alpha_deg)
            beta_rad = math.radians(beta_deg)
            gamma_rad = math.radians(gamma_deg)

            cos_a = math.cos(alpha_rad)
            cos_b = math.cos(beta_rad)
            cos_g = math.cos(gamma_rad)
            
            term = 1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g
            U_calc = a * b * c * math.sqrt(term)
            
            print(f"Equation: U = a*b*c * sqrt(1 - cos²α - cos²β - cos²γ + 2*cosα*cosβ*cosγ)")
            print(f"U = {a} * {b} * {c} * sqrt(1 - cos²({alpha_deg}) - cos²({beta_deg}) - cos²({gamma_deg}) + 2*cos({alpha_deg})*cos({beta_deg})*cos({gamma_deg}))")
            
        elif system == "monoclinic":
            a, b, c = data["a"], data["b"], data["c"]
            beta_deg = data["beta"]
            beta_rad = math.radians(beta_deg)
            
            U_calc = a * b * c * math.sin(beta_rad)
            
            print(f"Equation: U = a*b*c * sin(β)")
            print(f"U = {a} * {b} * {c} * sin({beta_deg})")

        elif system == "orthorhombic":
            a, b, c = data["a"], data["b"], data["c"]
            U_calc = a * b * c
            
            print(f"Equation: U = a * b * c")
            print(f"U = {a} * {b} * {c}")

        diff = abs(U_calc - U_reported)
        
        print(f"Calculated Volume = {U_calc:.2f} Å³")
        print(f"Reported Volume   = {U_reported:.2f} Å³")
        print(f"Difference        = {diff:.2f} Å³")
        
        if diff > 0.5: # A threshold to identify significant discrepancy
             print("Result: The reported volume is INCONSISTENT with the lattice parameters.\n")
        else:
             print("Result: The reported volume is consistent with the lattice parameters.\n")

    print("Conclusion: Dataset A shows a significant discrepancy between its lattice parameters and the reported unit cell volume.")


if __name__ == "__main__":
    check_crystal_data()