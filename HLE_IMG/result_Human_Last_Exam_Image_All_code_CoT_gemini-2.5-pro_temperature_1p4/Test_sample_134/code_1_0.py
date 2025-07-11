def name_the_molecule():
    """
    Identifies the molecule from its structure, calculates its chemical formula,
    and prints the information.
    """
    # 1. Identify and count the structural units from the image
    num_phenylene = 8
    num_ethynylene = 4
    num_butadiynylene = 4

    # 2. Define atoms per unit
    c_in_phenylene, h_in_phenylene = 6, 4
    c_in_ethynylene, h_in_ethynylene = 2, 0
    c_in_butadiynylene, h_in_butadiynylene = 4, 0

    # 3. Calculate total carbon atoms
    c_from_phenylene = num_phenylene * c_in_phenylene
    c_from_ethynylene = num_ethynylene * c_in_ethynylene
    c_from_butadiynylene = num_butadiynylene * c_in_butadiynylene
    total_c = c_from_phenylene + c_from_ethynylene + c_from_butadiynylene

    # 4. Calculate total hydrogen atoms
    total_h = num_phenylene * h_in_phenylene

    # 5. Print the name and analysis
    molecule_name = "Butadiyne-expanded cycloparaphenyleneethynylene"
    print(f"Molecule Name: {molecule_name}")
    print("\nThis molecule is a nanohoop macrocycle composed of alternating units.")

    # 6. Print the chemical formula calculation (the "equation")
    print("\nChemical Formula Calculation:")
    print("Carbon (C) atoms:")
    print(f"  From {num_phenylene} phenylene units: {num_phenylene} * {c_in_phenylene} = {c_from_phenylene}")
    print(f"  From {num_ethynylene} ethynylene units: {num_ethynylene} * {c_in_ethynylene} = {c_from_ethynylene}")
    print(f"  From {num_butadiynylene} butadiynylene units: {num_butadiynylene} * {c_in_butadiynylene} = {c_from_butadiynylene}")
    print(f"Total Carbon atoms = {c_from_phenylene} + {c_from_ethynylene} + {c_from_butadiynylene} = {total_c}")

    print("\nHydrogen (H) atoms:")
    print(f"  From {num_phenylene} phenylene units: {num_phenylene} * {h_in_phenylene} = {total_h}")
    print(f"Total Hydrogen atoms = {total_h}")

    print(f"\nFinal Chemical Formula: C{total_c}H{total_h}")

if __name__ == "__main__":
    name_the_molecule()