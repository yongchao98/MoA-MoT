def identify_molecule():
    """
    This function identifies the molecule in the image by analyzing its components,
    calculating its chemical formula, and providing its name.
    """
    # Number of each repeating unit in the macrocycle
    n_naphthylene = 6
    n_ethenylene = 6
    n_ethynylene = 12

    # Number of atoms in each unit
    c_naphthylene, h_naphthylene = 10, 6
    c_ethenylene, h_ethenylene = 2, 2
    c_ethynylene, h_ethynylene = 2, 0

    # Calculate total carbon and hydrogen atoms
    total_c = (n_naphthylene * c_naphthylene) + (n_ethenylene * c_ethenylene) + (n_ethynylene * c_ethynylene)
    total_h = (n_naphthylene * h_naphthylene) + (n_ethenylene * h_ethenylene) + (n_ethynylene * h_ethynylene)

    # Print the analysis and result
    print("Molecule Identification:")
    print("-" * 25)
    print("The molecule is composed of the following units:")
    print(f"- {n_naphthylene} units of 2,6-naphthalenediyl (C{c_naphthylene}H{h_naphthylene})")
    print(f"- {n_ethenylene} units of (Z)-ethenylene (C{c_ethenylene}H{h_ethenylene})")
    print(f"- {n_ethynylene} units of ethynylene (C{c_ethynylene})")
    print("\nChemical Formula Calculation:")
    print(f"Total Carbon = ({n_naphthylene} * {c_naphthylene}) + ({n_ethenylene} * {c_ethenylene}) + ({n_ethynylene} * {c_ethynylene}) = {n_naphthylene*c_naphthylene} + {n_ethenylene*c_ethenylene} + {n_ethynylene*c_ethynylene} = {total_c}")
    print(f"Total Hydrogen = ({n_naphthylene} * {h_naphthylene}) + ({n_ethenylene} * {h_ethenylene}) + ({n_ethynylene} * {h_ethynylene}) = {n_naphthylene*h_naphthylene} + {n_ethenylene*h_ethenylene} + {n_ethynylene*h_ethynylene} = {total_h}")
    print(f"\nFinal Chemical Formula: C{total_c}H{total_h}")

    print("\nName of the molecule:")
    print("A common name for this molecule is the C\u2089\u2086 nanohoop.")
    print("A descriptive name is: Cyclo-hexa[((2,6-naphthylene)diethynylene)-alt-((Z)-ethenylene)]")


if __name__ == '__main__':
    identify_molecule()