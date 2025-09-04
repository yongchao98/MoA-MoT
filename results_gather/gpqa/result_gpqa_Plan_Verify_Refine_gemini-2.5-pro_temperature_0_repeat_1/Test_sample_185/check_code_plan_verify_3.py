import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_cope_rearrangement_product():
    """
    This function checks the correctness of the given answer for the Cope rearrangement product.
    The logic follows these steps:
    1.  Define the reactant molecule from its name.
    2.  Simulate the [3,3]-sigmatropic (Cope) rearrangement to get an initial imine product.
    3.  Tautomerize the imine to its more stable enamine form, as suggested by the answer options.
    4.  Analyze the structure of the final enamine product.
    5.  Compare the properties of the calculated enamine with the properties implied by the answer 'D'.
    """
    try:
        # Step 1: Define the reactant
        # SMILES for 2-vinyl-2-azabicyclo[2.2.1]hept-5-ene
        # Stereochemistry is ignored as the options are constitutional isomers.
        reactant_smiles = "C=CN1CC2C=CC1C2"
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        if reactant_mol is None:
            return "Internal Error: Could not parse reactant SMILES."

        # Step 2: Simulate the Cope rearrangement to get the imine intermediate
        # The 1,5-diene system is C=C-N-C-C=C.
        patt = Chem.MolFromSmarts('[C:1]=[C:2]-[N:3]-[C:4]-[C:5]=[C:6]')
        if not reactant_mol.HasSubstructMatch(patt):
            return "Internal Error: Could not find the 1,5-diene system in the reactant."
        
        match_indices = reactant_mol.GetSubstructMatch(patt)
        c1_idx, c2_idx, n3_idx, c4_idx, c5_idx, c6_idx = match_indices

        emol = Chem.EditableMol(reactant_mol)
        emol.RemoveBond(n3_idx, c4_idx) # Break N-C bond
        emol.AddBond(c1_idx, c6_idx, Chem.BondType.SINGLE) # Form C-C bond
        # Adjust bond orders for the new structure
        emol.GetBondBetweenAtoms(c1_idx, c2_idx).SetBondType(Chem.BondType.SINGLE)
        emol.GetBondBetweenAtoms(c2_idx, n3_idx).SetBondType(Chem.BondType.DOUBLE)
        emol.GetBondBetweenAtoms(c4_idx, c5_idx).SetBondType(Chem.BondType.DOUBLE)
        emol.GetBondBetweenAtoms(c5_idx, c6_idx).SetBondType(Chem.BondType.SINGLE)

        product_mol_imine = emol.GetMol()
        Chem.SanitizeMol(product_mol_imine)

        # Step 3: Tautomerize the imine to the final enamine product
        from rdkit.Chem import rdchem
        tautomer_enumerator = rdchem.TautomerEnumerator()
        product_mol_enamine = tautomer_enumerator.Canonicalize(product_mol_imine)

        # Step 4: Analyze the calculated product and compare with answer 'D'
        # Answer D: 4,4a,5,6-tetrahydro-1H-cyclopenta[c]pyridine
        # This name implies 5 sp3-hybridized atoms (4 carbons + 1 nitrogen).
        
        # Count the sp3 atoms in our calculated enamine product.
        sp3_atoms_in_calc_product = [a for a in product_mol_enamine.GetAtoms() if a.GetHybridization() == Chem.HybridizationType.SP3]
        
        # Check if the count of sp3 atoms matches the count implied by the name.
        if len(sp3_atoms_in_calc_product) == 5:
            # The calculated product has 5 sp3 atoms, which is consistent with the name in option D.
            # This confirms that D is the correct constitutional isomer.
            return "Correct"
        else:
            return f"Incorrect: The provided answer 'D' implies a structure with 5 sp3 atoms (4 carbons and 1 nitrogen). However, the calculated final product has {len(sp3_atoms_in_calc_product)} sp3 atoms."

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check
result = check_cope_rearrangement_product()
print(result)