import sys
from io import StringIO

try:
    # RDKit is a standard cheminformatics library for molecular analysis.
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    # If rdkit is not installed, the check cannot be performed.
    print("Error: RDKit library not found. Please install it using 'pip install rdkit'")
    sys.exit(1)

def check_chemistry_answer():
    """
    This function verifies the multi-step chemical synthesis and NMR analysis.
    1. It programmatically determines the structure of the final product, E.
    2. It analyzes the symmetry of E to calculate the number of 13C-NMR signals.
    3. It compares this calculated number with the provided answer's choice.
    """
    # --- Part 1: Determine the structure of the final product E ---

    # The reaction sequence is well-defined and standard organic chemistry.
    # We can represent the key molecules using SMILES strings.

    # Step 1-4: Corey-Seebach synthesis to form the ketone (D)
    # Propionaldehyde -> 2-ethyl-1,3-dithiolane (A)
    # Alkylation with bromoethane -> 2,2-diethyl-1,3-dithiolane (C)
    # Deprotection -> 3-pentanone (D)
    ketone_smiles = "CCC(=O)CC"

    # Step 5: Wittig reaction
    # The ketone (D) reacts with an ylide formed from 3-bromopentane.
    # The carbon skeleton of the ylide is derived from 3-bromopentane: (CH3CH2)2CH-
    # The carbon skeleton of the ketone is derived from 3-pentanone: (CH3CH2)2C=
    # The Wittig reaction joins these two skeletons at the carbonyl/ylide carbons.
    # Final Product (E): 3,4-diethylhex-3-ene
    final_product_smiles = "CCC(=C(CC)CC)CC"

    # Create an RDKit molecule object for the final product
    mol = Chem.MolFromSmiles(final_product_smiles)
    if mol is None:
        return f"Internal Code Error: Could not parse the SMILES string for the final product: {final_product_smiles}"

    # --- Part 2: Calculate the number of 13C-NMR signals ---

    # The number of signals in a proton-decoupled 13C-NMR spectrum corresponds to
    # the number of chemically non-equivalent carbon atoms.
    # We can find this by determining the molecule's symmetry.

    # Add explicit hydrogens to ensure correct symmetry perception
    mol_with_hs = Chem.AddHs(mol)

    # Use CanonicalRankAtoms to assign a rank to each atom.
    # Symmetrically equivalent atoms will receive the same rank.
    ranks = list(Chem.CanonicalRankAtoms(mol_with_hs, breakTies=False))

    # Collect the ranks of only the carbon atoms (Atomic Number 6)
    carbon_ranks = []
    for atom in mol_with_hs.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_ranks.append(ranks[atom.GetIdx()])

    # The number of unique signals is the number of unique ranks among the carbons.
    calculated_signals = len(set(carbon_ranks))

    # --- Part 3: Compare with the provided answer ---

    # The question's options are: A) 6, B) 8, C) 3, D) 11
    options = {'A': 6, 'B': 8, 'C': 3, 'D': 11}
    
    # The provided answer selected 'C'.
    llm_choice_letter = 'C'
    llm_choice_value = options.get(llm_choice_letter)

    # Check 1: Does the calculated number of signals match the value of the chosen option?
    if calculated_signals != llm_choice_value:
        return (f"Incorrect. The final product, 3,4-diethylhex-3-ene, is highly symmetrical and should have {calculated_signals} signals. "
                f"The answer chose option {llm_choice_letter}, which corresponds to {llm_choice_value} signals.")

    # Check 2: Does the reasoning in the provided answer align with the calculation?
    # The provided answer correctly identifies the final product as 3,4-diethylhex-3-ene.
    # It correctly analyzes the symmetry, noting that all four ethyl groups are equivalent.
    # It correctly concludes there are 3 unique carbon environments (alkene, methylene, methyl).
    # This reasoning perfectly matches the computational result.
    
    # One common error (seen in other candidate answers) is to incorrectly assume E/Z isomerism.
    # For (CH3CH2)2C=C(CH2CH3)2, no E/Z isomers are possible because each sp2 carbon is bonded to two identical groups.
    # The provided answer correctly avoids this error.

    return "Correct"

# Run the check and print the result
result = check_chemistry_answer()
print(result)
