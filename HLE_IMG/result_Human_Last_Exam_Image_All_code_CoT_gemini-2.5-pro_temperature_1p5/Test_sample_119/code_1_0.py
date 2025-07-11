def solve_iupac_name():
    """
    This function analyzes the provided spectral data (MS, IR, 1H NMR, 13C NMR, DEPT, HSQC)
    to determine the structure and IUPAC name of the unknown organic compound.
    The analysis logic is summarized in the print statements.
    """

    # --- Step 1: Fragment Identification ---
    # From MS, IR, and NMR data, the key fragments are identified:
    # - A phenyl group (C6H5-)
    # - A primary amine group (-NH2)
    # - A methylene group (-CH2-)
    # - A methine group (-CH-)
    # - A methyl group (-CH3)
    
    # --- Step 2: Structure Assembly and Naming ---
    # Assembling these fragments consistent with all spectral data (especially NMR coupling)
    # leads to the structure: C6H5-CH2-CH(NH2)-CH3.
    # IUPAC Naming logic:
    # - The parent chain with the amine is 3 carbons long: propane.
    # - The amine is on position 2: propan-2-amine.
    # - The phenyl substituent is on position 1: 1-phenyl.
    
    print("Based on the comprehensive analysis of all provided spectra,")
    print("the determined molecular structure is C6H5-CH2-CH(NH2)-CH3.")
    
    name_part_1 = "1"
    name_part_2 = "-phenylpropan-"
    name_part_3 = "2"
    name_part_4 = "-amine"
    
    final_name = name_part_1 + name_part_2 + name_part_3 + name_part_4

    print("\nThe IUPAC name is constructed as follows:")
    print("Parent Chain: propane")
    print(f"Amine Group at Position: {name_part_3}")
    print(f"Phenyl Group at Position: {name_part_1}")
    
    print("\n--- Final IUPAC Name ---")
    print(f"The name is: {name_part_1}{name_part_2}{name_part_3}{name_part_4}")


# Execute the analysis to display the result.
solve_iupac_name()