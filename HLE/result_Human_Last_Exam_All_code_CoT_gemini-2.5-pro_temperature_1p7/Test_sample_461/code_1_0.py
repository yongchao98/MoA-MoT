def calculate_molar_mass():
    """
    Calculates the molar mass of the polycyclic hydrocarbon Dinosaurane (C29H52).
    """
    # Define atomic masses (g/mol)
    atomic_mass_carbon = 12.011
    atomic_mass_hydrogen = 1.008

    # Define the number of atoms in Dinosaurane (C29H52)
    num_carbon = 29
    num_hydrogen = 52

    # Calculate the molar mass
    molar_mass = (num_carbon * atomic_mass_carbon) + (num_hydrogen * atomic_mass_hydrogen)

    # Print the explanation and the full calculation equation
    print("The chemical formula for Dinosaurane is C29H52.")
    print("To find its molar mass, we sum the masses of its constituent atoms.")
    print(f"Molar Mass = (Carbon Atoms * Carbon Mass) + (Hydrogen Atoms * Hydrogen Mass)")
    # The final print statement includes each number used in the calculation, as requested.
    print(f"Molar Mass = ({num_carbon} * {atomic_mass_carbon}) + ({num_hydrogen} * {atomic_mass_hydrogen}) = {molar_mass:.4f} g/mol")

calculate_molar_mass()