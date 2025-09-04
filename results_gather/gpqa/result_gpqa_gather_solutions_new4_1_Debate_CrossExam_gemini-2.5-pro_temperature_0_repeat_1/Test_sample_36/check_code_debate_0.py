import sys

def check_nmr_signals():
    """
    This function verifies the number of 13C-NMR signals for the final product E.
    It uses the RDKit library to determine the structure of E and analyze its symmetry.
    """
    try:
        from rdkit import Chem
    except ImportError:
        # Instructions for the user if RDKit is not installed.
        return ("The 'rdkit' library is required to run this check. "
                "Please install it using 'pip install rdkit-pypi'.")

    # Step 1: Determine the structure of the final product, E.
    # The reaction sequence is:
    # 1. Propionaldehyde + EDT / BF3 -> A (2-ethyl-1,3-dithiolane)
    # 2. A + BuLi -> B (Lithiated carbanion)
    # 3. B + Bromoethane -> C (2,2-diethyl-1,3-dithiolane)
    # 4. C + HgCl2 / H2O / H+ -> D (3-pentanone)
    # 5. D + PPh3 / 3-bromopentane / BuLi -> E (Wittig reaction)
    # The Wittig reaction combines the carbon skeleton of 3-pentanone ((CH3CH2)2C=)
    # with the carbon skeleton of the ylide from 3-bromopentane (=C(CH2CH3)2).
    # The final product E is 3,4-diethylhex-3-ene.
    
    # The SMILES (Simplified Molecular Input Line Entry System) string for 3,4-diethylhex-3-ene is:
    smiles_E = "CCC(CC)=C(CC)CC"

    # Step 2: Create a molecule object from the SMILES string.
    mol = Chem.MolFromSmiles(smiles_E)
    if mol is None:
        return f"Error: Could not create a molecule from the SMILES string '{smiles_E}'."

    # Add explicit hydrogens to the molecule for a more accurate symmetry perception.
    mol_with_hs = Chem.AddHs(mol)

    # Step 3: Calculate the number of unique carbon environments.
    # RDKit's CanonicalRankAtoms function assigns a rank to each atom.
    # Atoms that are symmetrically equivalent will have the same rank.
    # We count the number of unique ranks for carbon atoms only.
    ranks = Chem.CanonicalRankAtoms(mol_with_hs, breakTies=False)
    
    carbon_ranks = set()
    for atom in mol_with_hs.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Atomic number for Carbon is 6
            carbon_ranks.add(ranks[atom.GetIdx()])
            
    calculated_signals = len(carbon_ranks)

    # Step 4: Compare the calculated result with the provided answer.
    # The question's correct option is A, which corresponds to 3 signals.
    expected_signals = 3
    
    if calculated_signals == expected_signals:
        return "Correct"
    else:
        return (f"Incorrect. The final product is 3,4-diethylhex-3-ene. "
                f"A computational analysis of its structure reveals {calculated_signals} unique carbon environments, "
                f"which corresponds to {calculated_signals} 13C-NMR signals. "
                f"The provided answer states there are {expected_signals} signals.")

# Execute the check and print the result.
# Note: This block requires the 'rdkit-pypi' package to be installed.
try:
    result = check_nmr_signals()
    print(result)
except Exception as e:
    print(f"An error occurred during execution: {e}")
