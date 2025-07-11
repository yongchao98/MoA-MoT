def calculate_rarest_noble_gas_percentage():
    """
    Calculates and prints the abundance of the rarest noble gas (Radon)
    as a percentage of all terrestrial matter.
    """
    gas_name = "Radon"
    
    # Estimated total mass of Radon in the entire Earth's atmosphere at any given time, in kg.
    # Sources estimate this to be on the order of a few thousand kilograms. We'll use 2000 kg.
    mass_of_gas_kg = 2000.0
    
    # Mass of the Earth in kg.
    mass_of_earth_kg = 5.972e24
    
    # Calculate the percentage.
    percentage = (mass_of_gas_kg / mass_of_earth_kg) * 100
    
    print(f"The rarest noble gas on Earth is {gas_name}.")
    print("Its abundance as a percentage of all terrestrial matter is extremely low.")
    print("The calculation is based on the estimated total mass of the gas on Earth divided by the total mass of the Earth.")
    print("\nCalculation:")
    # The prompt requires printing each number in the final equation
    print(f"({mass_of_gas_kg} kg / {mass_of_earth_kg} kg) * 100 = {percentage:.4g}%")

calculate_rarest_noble_gas_percentage()