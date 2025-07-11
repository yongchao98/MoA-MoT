def calculate_rarest_noble_gas_percentage():
    """
    Calculates the percentage of the rarest noble gas (Radon)
    in all terrestrial matter.
    """
    # Scientific data for the calculation. These values are stable and
    # representative for the year 2002.
    gas_name = "Radon (Rn)"
    
    # Mass of the Earth in kilograms
    mass_earth_kg = 5.972e24
    
    # Mass of the Earth's atmosphere in kilograms
    mass_atmosphere_kg = 5.15e18
    
    # Mass fraction of Radon in the atmosphere. This is an extremely small number.
    radon_mass_fraction_in_atmosphere = 6e-17

    print(f"The rarest noble gas on Earth is {gas_name}.")
    print("This script calculates its approximate percentage of all terrestrial matter.")
    print("-" * 30)

    # Step 1: Calculate the total mass of Radon in the atmosphere
    total_mass_radon_kg = radon_mass_fraction_in_atmosphere * mass_atmosphere_kg
    
    # Step 2: Calculate the final percentage relative to Earth's total mass
    percentage_of_earth = (total_mass_radon_kg / mass_earth_kg) * 100

    # Print the equation and the final result
    print("Equation: ((Radon Mass Fraction in Atmosphere * Mass of Atmosphere) / Mass of Earth) * 100")
    print(f"Percentage = (({radon_mass_fraction_in_atmosphere} * {mass_atmosphere_kg}) / {mass_earth_kg}) * 100")
    print(f"Final Percentage: {percentage_of_earth:.4e} %")

calculate_rarest_noble_gas_percentage()

<<<Radon (Rn)>>>