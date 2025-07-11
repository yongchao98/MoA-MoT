def name_molecule():
    """
    This function provides the name and molecular formula for the given molecule.
    """
    molecule_name = "Cyclo[4](1,4-phenylene)-(E)-hex-3-en-1,5-diyne"
    
    # Details of the molecular formula calculation
    c_phenylene = 4 * 6
    c_linker = 4 * 6
    total_c = c_phenylene + c_linker
    
    h_phenylene = 4 * 4
    h_linker = 4 * 2
    total_h = h_phenylene + h_linker
    
    print(f"Molecule Name: {molecule_name}")
    print("\nDerivation of Molecular Formula (C48H24):")
    print("Carbon (C) atoms:")
    print(f"From 4 phenylene rings: 4 * 6 = {c_phenylene}")
    print(f"From 4 linker chains: 4 * 6 = {c_linker}")
    print(f"Total Carbon atoms: {c_phenylene} + {c_linker} = {total_c}")
    
    print("\nHydrogen (H) atoms:")
    print(f"From 4 phenylene rings: 4 * 4 = {h_phenylene}")
    print(f"From 4 linker chains: 4 * 2 = {h_linker}")
    print(f"Total Hydrogen atoms: {h_phenylene} + {h_linker} = {total_h}")
    
    print(f"\nFinal Molecular Formula: C{total_c}H{total_h}")

name_molecule()