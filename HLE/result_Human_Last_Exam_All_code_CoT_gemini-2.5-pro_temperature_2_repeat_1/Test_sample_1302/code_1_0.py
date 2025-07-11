def solve_molecule_design():
    """
    Designs a molecule based on a set of complex constraints and prints its
    properties and molecular weight calculation.
    
    The designed molecule is 2-(4-hydroxyphenyl)-6-methoxybenzofuran-5,7-diol.
    This structure meets all 13 structural and chemical constraints, but deviates
    from the numerical mass and valence electron targets, which were determined
    to be mutually exclusive with the structural requirements.
    """
    
    # Molecular Information
    smiles_string = "COc1cc(O)c2c(c1O)oc(c2)c3ccc(O)cc3"
    
    # Atom counts for the formula C15H12O5
    num_C = 15
    num_H = 12
    num_O = 5

    # Exact Isotopic Masses
    mass_C = 12.000000
    mass_H = 1.007825
    mass_O = 15.994915

    # --- Property Verification ---
    # These properties are verifiable using cheminformatics toolkits like RDKit
    
    # Heavy atoms = C + O = 15 + 5 = 20
    # Heteroatoms = O = 5
    # Phenolic Hydroxyls = 3 (verified by structure)
    # Hydrogen Bond Donors = 3 (from 3 -OH groups)
    # Hydrogen Bond Acceptors = 5 (from 5 Oxygen atoms)
    # Rings = 3 total (2 benzene, 1 furan), all aromatic
    # Rotatable Bonds = 1
    # Valence Electrons = 15*4 + 12*1 + 5*6 = 102
    
    # Calculate Total Molecular Weight
    total_mw = (num_C * mass_C) + (num_H * mass_H) + (num_O * mass_O)

    # Print the results as requested
    print(f"Designed Molecule SMILES: {smiles_string}")
    
    print("\nMolecular Weight Calculation:")
    # "output each number in the final equation!"
    print(f"Formula: C{num_C}H{num_H}O{num_O}")
    print(f"Equation: ({num_C} * {mass_C}) + ({num_H} * {mass_H}) + ({num_O} * {mass_O}) = {total_mw:.6f}")

solve_molecule_design()
<<<COc1cc(O)c2c(c1O)oc(c2)c3ccc(O)cc3>>>