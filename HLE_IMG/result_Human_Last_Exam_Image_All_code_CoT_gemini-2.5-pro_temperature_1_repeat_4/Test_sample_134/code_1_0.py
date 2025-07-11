def identify_molecule():
    """
    This function identifies a specific carbon nanobelt by calculating its
    chemical formula from its constituent parts and then provides its name.
    """
    # Step 1 & 2: Define the number of each structural unit identified from the image.
    num_phenylene = 12
    num_ethenylene = 6
    num_ethynylene = 8

    # Define the number of C and H atoms per unit.
    atoms_in_phenylene = {'C': 6, 'H': 4}
    atoms_in_ethenylene = {'C': 2, 'H': 2}
    atoms_in_ethynylene = {'C': 2, 'H': 0}

    # Step 3: Calculate the total number of carbon and hydrogen atoms.
    total_C = (num_phenylene * atoms_in_phenylene['C']) + \
              (num_ethenylene * atoms_in_ethenylene['C']) + \
              (num_ethynylene * atoms_in_ethynylene['C'])
              
    total_H = (num_phenylene * atoms_in_phenylene['H']) + \
              (num_ethenylene * atoms_in_ethenylene['H']) + \
              (num_ethynylene * atoms_in_ethynylene['H'])

    # Print the breakdown of the calculation for transparency.
    print("Calculation of the Chemical Formula:")
    print(f"Carbon atoms = ({num_phenylene} * {atoms_in_phenylene['C']}) + ({num_ethenylene} * {atoms_in_ethenylene['C']}) + ({num_ethynylene} * {atoms_in_ethynylene['C']}) = {num_phenylene * atoms_in_phenylene['C']} + {num_ethenylene * atoms_in_ethenylene['C']} + {num_ethynylene * atoms_in_ethynylene['C']} = {total_C}")
    print(f"Hydrogen atoms = ({num_phenylene} * {atoms_in_phenylene['H']}) + ({num_ethenylene} * {atoms_in_ethenylene['H']}) = {num_phenylene * atoms_in_phenylene['H']} + {num_ethenylene * atoms_in_ethenylene['H']} = {total_H}")
    print("\n------------------------------------\n")

    # Step 4: Print the final name of the molecule.
    print(f"Molecule Name: C{total_C}H{total_H} carbon nanobelt")

if __name__ == "__main__":
    identify_molecule()