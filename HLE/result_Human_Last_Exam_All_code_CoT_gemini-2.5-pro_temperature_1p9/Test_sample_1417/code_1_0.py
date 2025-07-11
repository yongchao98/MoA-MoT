import operator

def solve():
    """
    Calculates the optimum material for a particle detector's cooling tubes based
    on the unique parameter of maximizing radiation length.
    """
    # The unique parameter for a particle detector is minimizing interaction with
    # particles, which means maximizing the radiation length (X0).
    # For pipes of similar dimensions, we need to compare the radiation length
    # in cm, calculated as: X0_cm = X0_g_cm2 / density_g_cm3
    # Data is from the Particle Data Group (PDG).
    materials = {
        'Titanium': {'density_g_cm3': 4.51, 'X0_g_cm2': 16.0},
        'Aluminium': {'density_g_cm3': 2.70, 'X0_g_cm2': 24.01},
        '316 Stainless Steel': {'density_g_cm3': 8.00, 'X0_g_cm2': 13.9}, # Approximated with Iron values
        'Copper': {'density_g_cm3': 8.96, 'X0_g_cm2': 12.86},
        'Nickel': {'density_g_cm3': 8.91, 'X0_g_cm2': 12.49}
    }

    results = {}

    print("The unique parameter is radiation length (X0). A longer X0 is better.")
    print("Calculating X0 in cm (X0_g_cm2 / density) for each material:")
    print("---------------------------------------------------------------")

    for name, properties in materials.items():
        density = properties['density_g_cm3']
        x0_mass_thickness = properties['X0_g_cm2']
        # The equation for the final calculation
        x0_length = x0_mass_thickness / density
        results[name] = x0_length
        # Outputting each number in the final equation as requested
        print(f"{name}: {x0_mass_thickness} / {density} = {x0_length:.2f} cm")

    # Find the material with the maximum radiation length in cm.
    optimum_choice, max_x0 = max(results.items(), key=operator.itemgetter(1))

    print("---------------------------------------------------------------")
    print(f"The optimum material is the one with the longest radiation length.")
    print(f"Winner: {optimum_choice} with X0 = {max_x0:.2f} cm")

solve()
<<<B>>>