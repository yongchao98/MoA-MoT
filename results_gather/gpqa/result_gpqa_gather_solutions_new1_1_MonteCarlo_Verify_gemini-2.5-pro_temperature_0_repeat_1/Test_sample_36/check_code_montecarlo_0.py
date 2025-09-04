import sys
try:
    from rdkit import Chem
except ImportError:
    # If rdkit is not installed, the check cannot be performed.
    # This is a common library in cheminformatics.
    # To install: pip install rdkit-pypi
    print("RDKit library not found. Please install it to run this check.")
    sys.exit(1)

def check_nmr_signals():
    """
    This function verifies the number of 13C-NMR signals for the final product of the given reaction sequence.
    """
    # Step 1: Define the final product based on the reaction analysis.
    # The reaction sequence is a Corey-Seebach synthesis followed by a Wittig reaction.
    # Reactant D (ketone): 3-pentanone -> CCC(=O)CC
    # Wittig reagent from 3-bromopentane -> ylide is (CH3CH2)2C=PPh3
    # The final product E is 3,4-diethylhex-3-ene.
    # Its structure is (CH3CH2)2C=C(CH2CH3)2.
    smiles_E = "CCC(CC)=C(CC)CC"

    # Step 2: Create an RDKit molecule object from the SMILES string.
    mol = Chem.MolFromSmiles(smiles_E)
    if mol is None:
        return "Error: Could not parse the SMILES string for the final product."

    # A critical point is that no E/Z isomers are possible for this molecule,
    # as each carbon of the double bond is attached to two identical groups (ethyl groups).
    # Therefore, analyses considering E/Z isomers are incorrect.

    # Step 3: Calculate the number of chemically non-equivalent carbons.
    # RDKit's canonical atom ranking algorithm identifies symmetrically equivalent atoms.
    # Atoms with the same rank are considered equivalent.
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=True))

    # Create a set to store the unique ranks of only the carbon atoms.
    carbon_ranks = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Atomic number for Carbon is 6
            rank_for_this_atom = ranks[atom.GetIdx()]
            carbon_ranks.add(rank_for_this_atom)

    calculated_signals = len(carbon_ranks)

    # Step 4: Compare the calculated number of signals with the provided answer.
    # The question's options are: A) 8, B) 3, C) 11, D) 6
    options = {'A': 8, 'B': 3, 'C': 11, 'D': 6}
    
    # The provided answer is <<<B>>>, which corresponds to a value of 3.
    provided_answer_choice = 'B'
    provided_answer_value = options[provided_answer_choice]

    # Step 5: Return the result of the check.
    if calculated_signals == provided_answer_value:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer states there are {provided_answer_value} signals (Option {provided_answer_choice}). "
                f"However, a computational analysis of the final product, 3,4-diethylhex-3-ene, shows it has {calculated_signals} "
                f"unique carbon environments. The correct number of signals should be {calculated_signals}.")

# Run the check and print the result.
result = check_nmr_signals()
print(result)