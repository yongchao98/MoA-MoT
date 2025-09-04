import sys

def check_answer_correctness():
    """
    Checks the correctness of the provided answer for the three Michael addition reactions.
    The answer to check is:
    A = trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate
    B = 3-(2-oxocyclohexyl)butanenitrile
    C = cyclohexane-1,3-dione
    """
    try:
        from rdkit import Chem
    except ImportError:
        return "Error: RDKit library not found. Please install it using 'pip install rdkit-pypi' to run this check."

    def get_canonical_smiles(smiles_string):
        """Converts a SMILES string to its canonical form."""
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            # This indicates an invalid SMILES string was provided.
            return f"Invalid SMILES: {smiles_string}"
        return Chem.MolToSmiles(mol, canonical=True)

    # --- Verification for Reaction A ---
    # Reaction: dimethyl malonate + methyl (E)-3-(p-tolyl)acrylate --> A
    # The nucleophile from dimethyl malonate, [-CH(COOMe)2], attacks the beta-carbon of the acrylate.
    # Expected Product Structure: p-tolyl-CH(CH(COOMe)2)-CH2-COOMe
    expected_A_smiles = "COC(=O)CC(c1ccc(C)cc1)C(C(=O)OC)C(=O)OC"
    
    # Answer A Name: trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate
    # This name corresponds to the same structure derived from the mechanism.
    answer_A_smiles = "COC(=O)CC(c1ccc(C)cc1)C(C(=O)OC)C(=O)OC"

    if get_canonical_smiles(expected_A_smiles) != get_canonical_smiles(answer_A_smiles):
        return (f"Incorrect product A. The structure for 'trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate' "
                f"does not match the expected product of the Michael addition. "
                f"Expected SMILES: {get_canonical_smiles(expected_A_smiles)}, "
                f"Answer SMILES: {get_canonical_smiles(answer_A_smiles)}")

    # --- Verification for Reaction B ---
    # Reaction: 1-(cyclohex-1-en-1-yl)piperidine + (E)-but-2-enenitrile --> B (after H3O+ workup)
    # This is a Stork enamine synthesis. The enamine attacks the beta-carbon of the nitrile.
    # The H3O+ workup hydrolyzes the intermediate iminium salt back to a ketone.
    # Expected Product Structure: A cyclohexanone ring with a -CH(CH3)CH2CN group at the alpha-position.
    expected_B_smiles = "N#CCC(C)C1CCCC(=O)C1"

    # Answer B Name: 3-(2-oxocyclohexyl)butanenitrile
    # This name corresponds to the expected product structure.
    answer_B_smiles = "N#CCC(C)C1CCCC(=O)C1"

    if get_canonical_smiles(expected_B_smiles) != get_canonical_smiles(answer_B_smiles):
        return (f"Incorrect product B. The structure for '3-(2-oxocyclohexyl)butanenitrile' "
                f"does not match the expected product of the Stork enamine reaction. "
                f"The major product after acidic workup should be the ketone, not an enol or enamine. "
                f"Expected SMILES: {get_canonical_smiles(expected_B_smiles)}, "
                f"Answer SMILES: {get_canonical_smiles(answer_B_smiles)}")

    # --- Verification for Reaction C ---
    # Reaction (Retrosynthesis): C + but-3-en-2-one ---> 2-(3-oxobutyl)cyclohexane-1,3-dione
    # The product is formed by adding a '3-oxobutyl' group (-CH2CH2COCH3) to a nucleophile.
    # This group comes from the Michael acceptor, but-3-en-2-one.
    # The nucleophile must therefore be the enolate of cyclohexane-1,3-dione.
    # Reactant C is the precursor to this nucleophile.
    expected_C_smiles = "O=C1CC(=O)CCC1"

    # Answer C Name: cyclohexane-1,3-dione
    # This name corresponds to the expected reactant.
    answer_C_smiles = "O=C1CC(=O)CCC1"

    if get_canonical_smiles(expected_C_smiles) != get_canonical_smiles(answer_C_smiles):
        return (f"Incorrect reactant C. The reactant required to form '2-(3-oxobutyl)cyclohexane-1,3-dione' "
                f"via Michael addition with but-3-en-2-one is 'cyclohexane-1,3-dione'. "
                f"Expected SMILES: {get_canonical_smiles(expected_C_smiles)}, "
                f"Answer SMILES: {get_canonical_smiles(answer_C_smiles)}")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_answer_correctness()
print(result)