import operator

def solve_detector_material():
    """
    This function determines the optimal material for a particle detector's
    cooling system by maximizing the radiation length, a key parameter for
    minimizing particle interaction.
    """

    # Data for materials:
    # 'density' is in g/cm^3
    # 'rad_length_mass' is the radiation length in g/cm^2
    # Data is sourced from the Particle Data Group (PDG) where available.
    materials = {
        'Titanium': {'density': 4.51, 'rad_length_mass': 16.0},
        'Aluminium': {'density': 2.70, 'rad_length_mass': 24.01},
        '316 Stainless Steel': {'density': 8.00, 'rad_length_mass': 13.9}, # Approximation for Iron/Steel
        'Copper': {'density': 8.96, 'rad_length_mass': 12.86},
        'Nickle': {'density': 8.91, 'rad_length_mass': 12.49}
    }

    # This dictionary will store the calculated radiation length in cm for each material
    results_cm = {}

    print("The unique parameter to maximize is the radiation length (X₀) to minimize particle interaction.")
    print("Calculating X₀ in cm for each material using the formula: X₀ [cm] = X₀ [g/cm^2] / Density [g/cm^3]\n")

    for name, properties in materials.items():
        density = properties['density']
        rad_length_mass = properties['rad_length_mass']

        # Calculate the radiation length in units of length (cm)
        rad_length_cm = rad_length_mass / density
        results_cm[name] = rad_length_cm

        # Output the calculation for each material as requested
        print(f"For {name}:")
        print(f"  The equation is: {rad_length_mass} / {density} = {rad_length_cm:.2f}")
        print(f"  The radiation length is {rad_length_cm:.2f} cm.\n")


    # Find the material with the maximum radiation length in cm
    # The max function with a key operates on the dictionary's items (key-value pairs)
    # and operator.itemgetter(1) tells it to use the value (the second element) for comparison.
    optimum_choice, max_rad_length = max(results_cm.items(), key=operator.itemgetter(1))

    print("--- Conclusion ---")
    print(f"The material with the longest radiation length is '{optimum_choice}'.")
    print(f"This makes it the most 'transparent' to particles and therefore the optimum choice for the cooling system tubes.")

solve_detector_material()
<<<B>>>