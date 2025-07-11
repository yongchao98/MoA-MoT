def identify_molecule():
    """
    Identifies the molecule from the image and prints its name and formula.
    The molecule is a cyclic hexamer of para-phenylene and hexa-3-en-1,5-diynylene units.
    """
    # Define the components of the name
    cycle_size = 6
    phenylene_unit = "para-phenylene"
    linker_unit = "hexa-3-en-1,5-diynylene"
    
    # Construct the descriptive name
    molecule_name = f"Cyclo[{cycle_size}]{phenylene_unit}({linker_unit})"
    
    # Define the atom counts for the formula
    carbon_atoms = 72
    hydrogen_atoms = 36
    
    # Construct the molecular formula string
    molecular_formula = f"C{carbon_atoms}H{hydrogen_atoms}"

    print(f"Molecule Name: {molecule_name}")
    print(f"Molecular Formula: {molecular_formula}")

identify_molecule()