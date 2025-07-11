import sys

def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling system
    by maximizing the radiation length.
    """
    # Data for materials:
    # Density is in g/cm^3
    # Radiation_Length_Mass is the mass-normalized radiation length (X0) in g/cm^2
    materials = {
        "Titanium":         {"density": 4.51, "radiation_length_mass": 16.1},
        "Aluminium":        {"density": 2.70, "radiation_length_mass": 24.01},
        "316 Stainless Steel": {"density": 8.00, "radiation_length_mass": 13.84}, # Approximated with Iron's X0
        "Copper":           {"density": 8.96, "radiation_length_mass": 12.86},
        "Nickle":           {"density": 8.91, "radiation_length_mass": 12.93}
    }

    print("The unique parameter for a particle detector is to minimize particle interaction.")
    print("This is achieved by maximizing the radiation length (X0) of the material.\n")
    print("Calculating radiation length in cm for each material...")
    print("Formula: Radiation Length [cm] = Radiation Length [g/cm^2] / Density [g/cm^3]\n")

    optimal_material = None
    max_rad_length = -1

    for name, properties in materials.items():
        density = properties["density"]
        rad_length_mass = properties["radiation_length_mass"]
        
        # Calculate the radiation length in cm
        rad_length_cm = rad_length_mass / density
        
        # Output the calculation for each material
        print(f"For {name}:")
        print(f"  Radiation Length [cm] = {rad_length_mass} / {density} = {rad_length_cm:.2f} cm")
        
        if rad_length_cm > max_rad_length:
            max_rad_length = rad_length_cm
            optimal_material = name
            
    print("\n--------------------------------------------------")
    print(f"The material with the maximum radiation length is {optimal_material}.")
    print("Therefore, it is the optimum choice to minimize particle interaction.")
    print("--------------------------------------------------")

if __name__ == '__main__':
    find_optimal_material()