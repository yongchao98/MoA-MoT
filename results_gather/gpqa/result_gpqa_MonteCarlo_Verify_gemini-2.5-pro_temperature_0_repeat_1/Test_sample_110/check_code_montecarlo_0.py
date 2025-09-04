import sys

def check_chemistry_answer():
    """
    Checks the correctness of the proposed answer for two Michael addition reactions.
    It uses the RDKit library to compare chemical structures.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return "Error: RDKit library not found. Please install it (`pip install rdkit-pypi`) to run this verification code."

    def get_canonical_smiles(smiles_string):
        """Converts a SMILES string to its canonical form, ignoring stereochemistry for this problem."""
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            # This indicates an issue with the SMILES string itself.
            return f"InvalidSMILES_{smiles_string}"
        # The question does not specify stereochemistry, so we remove it for a robust comparison of connectivity.
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol, canonical=True)

    # --- Reaction A Verification ---
    # Reaction: 2-ethyl-2,6-dimethylcyclohexan-1-one + ethyl acrylate (t-BuOK)
    # Principle: The bulky base t-BuOK deprotonates the only available alpha-proton at the less-substituted C6 position.
    # The resulting enolate attacks the beta-carbon of ethyl acrylate in a Michael addition.
    
    # Structure of 2-ethyl-2,6-dimethylcyclohexan-1-one: CCC1(C)CCCC(C)C1=O
    # The product is formed by creating a bond between C6 and the beta-carbon of ethyl acrylate.
    # The methyl group at C6 remains.
    smiles_product_A_from_mechanism = "CCC1(C)CCCC(C(C)CCC(=O)OCC)C1=O"
    
    # Proposed product A name from option B: ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate
    # Let's build the SMILES from this name to verify its structure.
    # The name implies the cyclohexyl ring is a substituent on the propanoate chain.
    # Numbering the ring from the attachment point (C1'):
    # C1' has a methyl group.
    # C2' is the carbonyl (oxo).
    # C3' has an ethyl and a methyl group.
    # This matches the connectivity derived from the mechanism.
    smiles_product_A_from_name = "CCOC(=O)CCC1(C)C(=O)C(CC)(C)CCC1"

    canonical_A_mechanism = get_canonical_smiles(smiles_product_A_from_mechanism)
    canonical_A_name = get_canonical_smiles(smiles_product_A_from_name)

    if canonical_A_mechanism != canonical_A_name:
        return (f"Incorrect. The structure of product A is wrong.\n"
                f"Based on the mechanism (enolate formation at C6), the product should have the canonical SMILES: {canonical_A_mechanism}.\n"
                f"The product named in the answer ('ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate') has the canonical SMILES: {canonical_A_name}.\n"
                f"These structures do not match.")

    # --- Reaction B Verification ---
    # Reaction: 1-nitropropane + (E)-but-2-enenitrile (KOH)
    # Principle: KOH deprotonates the acidic alpha-carbon of 1-nitropropane.
    # The resulting nitronate anion attacks the beta-carbon of but-2-enenitrile.
    
    # Product formed by bonding alpha-carbon of nitropropane to beta-carbon of the nitrile.
    # Structure: CH3-CH2-CH(NO2)-CH(CH3)-CH2-CN
    smiles_product_B_from_mechanism = "CCC(C(C)CC#N)[N+](=O)[O-]"
    
    # Proposed product B name from option B: 3-methyl-4-nitrohexanenitrile
    # Let's build the SMILES from this name.
    # Main chain is hexanenitrile (6 carbons, CN is C1).
    # C3 has a methyl group.
    # C4 has a nitro group.
    # Structure: NC-CH2-CH(CH3)-CH(NO2)-CH2-CH3
    smiles_product_B_from_name = "N#CCC(C)C(C(C)[N+](=O)[O-])C" # A different way to write it, canonicalization will solve it.
    # A more direct translation is the same as the mechanism one:
    smiles_product_B_from_name_direct = "CCC(C(C)CC#N)[N+](=O)[O-]"

    canonical_B_mechanism = get_canonical_smiles(smiles_product_B_from_mechanism)
    canonical_B_name = get_canonical_smiles(smiles_product_B_from_name_direct)

    if canonical_B_mechanism != canonical_B_name:
        return (f"Incorrect. The structure of product B is wrong.\n"
                f"Based on the mechanism, the product should have the canonical SMILES: {canonical_B_mechanism}.\n"
                f"The product named in the answer ('3-methyl-4-nitrohexanenitrile') has the canonical SMILES: {canonical_B_name}.\n"
                f"These structures do not match.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_chemistry_answer()
print(result)
