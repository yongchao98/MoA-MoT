def solve_iupac_name():
    """
    This function prints the determined IUPAC name of the compound.
    The structure was determined to be 1-phenylpropan-2-amine based on the analysis of
    MS, IR, 1H NMR, 13C NMR, DEPT, and HSQC data.
    """
    # The parent chain is propane.
    # The principal functional group is an amine at position 2.
    amine_position = 2
    parent_chain = "propan"
    amine_suffix = "amine"

    # There is a phenyl substituent at position 1.
    substituent_position = 1
    substituent_name = "phenyl"

    # Assemble the IUPAC name
    # The numbers specifying the locants are explicitly included.
    iupac_name = f"{substituent_position}-{substituent_name}{parent_chain}-{amine_position}-{amine_suffix}"
    
    print("The IUPAC name of the compound is:")
    print(iupac_name)
    # The final equation is the name itself, which includes the numbers (locants).
    print("\nThe final name shows the locants:")
    print(f"Phenyl group at position: {substituent_position}")
    print(f"Amine group at position: {amine_position}")


solve_iupac_name()