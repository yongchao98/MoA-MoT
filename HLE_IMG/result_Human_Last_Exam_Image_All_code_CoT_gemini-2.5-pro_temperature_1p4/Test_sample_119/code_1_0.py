def solve_structure_and_name():
    """
    This function analyzes the provided spectroscopic data to determine the IUPAC name
    of the unknown compound and prints the analysis.
    """
    
    # Based on the step-by-step analysis of the provided spectra
    # MS, IR, 1H NMR, 13C NMR, DEPT, and HSQC.
    
    compound_name = "1-phenylpropan-2-amine"
    
    # The name is constructed based on IUPAC rules:
    # 1. Identify substituents and their positions.
    # 2. Identify the parent chain.
    # 3. Identify the principal functional group and its position.
    
    substituent_name = "phenyl"
    substituent_position = 1
    
    parent_chain_name = "propan"
    parent_chain_length = 3
    
    principal_group_name = "amine"
    principal_group_position = 2

    print(f"The determined IUPAC name of the compound is: {compound_name}")
    print("\nThis name is derived from the following components:")
    print(f"Substituent: {substituent_name} at position {substituent_position}")
    print(f"Parent Chain: {parent_chain_name} (length {parent_chain_length})")
    print(f"Principal Group: {principal_group_name} at position {principal_group_position}")
    
    print("\nThe numbers that form the IUPAC name locants are:")
    # Final 'equation' as requested by the prompt, showing the numbers in the name
    final_equation_numbers = f"{substituent_position}-{substituent_name}{parent_chain_name}-{principal_group_position}-{principal_group_name}"
    print(f"Position of phenyl group: {substituent_position}")
    print(f"Position of amine group: {principal_group_position}")
    print(f"Final name construction: {final_equation_numbers}")


solve_structure_and_name()
