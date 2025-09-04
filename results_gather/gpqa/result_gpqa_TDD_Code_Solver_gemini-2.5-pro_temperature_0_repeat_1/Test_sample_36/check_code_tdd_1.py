import sys

def check_correctness():
    """
    This function checks the correctness of a given answer to a multi-step organic chemistry synthesis problem.
    
    The problem is:
    Identify the number of 13C-NMR signals produced by the final product, E, resulting from the series of reactions shown below.
    Propionaldehyde + EDT / BF3 ---> A
    A + BuLi ---> B
    B + Bromoethane ---> C
    C + HgCl2 / H2O / H+ ---> D
    D + PPh3 / 3-bromopentane / BuLi ---> E

    The provided answer is:
    The final product is 3,4-diethylhex-3-ene, which has 3 unique carbon environments, corresponding to 3 signals (Option A).
    
    This checker will:
    1. Verify the identity of the final product E based on the reaction sequence.
    2. Calculate the number of 13C-NMR signals for product E using molecular symmetry.
    3. Compare the calculated result with the provided answer.
    """
    
    # First, check if RDKit is installed, as it is essential for this analysis.
    try:
        from rdkit import Chem
        from rdkit.Chem import rdmolfiles
    except ImportError:
        return "Execution Error: The RDKit library is required to perform this check. Please install it using 'pip install rdkit-pypi'."

    # --- Analysis of the Reaction Sequence ---

    # Step 1-4: Formation of the ketone (Product D)
    # Propionaldehyde (CCC=O) is protected as a dithiane (A: 2-ethyl-1,3-dithiane).
    # The dithiane is deprotonated and alkylated with bromoethane (C: 2,2-diethyl-1,3-dithiane).
    # The dithiane is deprotected to reveal the ketone (D: 3-pentanone or diethyl ketone).
    # SMILES for D: CCC(=O)CC
    
    # Step 5: Wittig Reaction to form the final alkene (Product E)
    # Ketone D (3-pentanone) reacts with a Wittig ylide.
    # The ylide is formed from 3-bromopentane (CCC(Br)CC) + PPh3, then BuLi.
    # The ylide is (CH3CH2)2C=PPh3.
    # The reaction is (Et)2C=O + (Et)2C=PPh3 -> (Et)2C=C(Et)2 + Ph3P=O.
    # The final product E is 3,4-diethylhex-3-ene.
    smiles_E = "CCC(=C(CC)CC)CC"
    
    # --- Verification of the Answer ---

    # Claim 1: The final product is 3,4-diethylhex-3-ene.
    # Our analysis confirms this structure, represented by the SMILES string `smiles_E`.
    
    # Claim 2: The product has 3 unique carbon environments (3 13C-NMR signals).
    # We will calculate this using RDKit's symmetry perception tools.

    def count_carbon_nmr_signals(smiles_string):
        """
        Calculates the number of unique 13C NMR signals for a molecule
        by finding the number of symmetrically non-equivalent carbon atoms.
        """
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return None # Error handling for invalid SMILES
            
        # Add explicit hydrogens to ensure the graph is complete for symmetry perception
        mol_with_h = Chem.AddHs(mol)
        
        # CanonicalRankAtoms assigns a rank to each atom. Symmetrically equivalent atoms get the same rank.
        ranks = list(rdmolfiles.CanonicalRankAtoms(mol_with_h, breakTies=False))
        
        # We are interested in the ranks of carbon atoms only.
        carbon_ranks = set()
        for i, atom in enumerate(mol_with_h.GetAtoms()):
            if atom.GetAtomicNum() == 6:  # It's a Carbon atom
                carbon_ranks.add(ranks[i])
        
        return len(carbon_ranks)

    calculated_signals = count_carbon_nmr_signals(smiles_E)
    
    if calculated_signals is None:
        return f"Internal Error: Could not process the SMILES string for the final product E: {smiles_E}"

    # The provided answer corresponds to option A, which is 3 signals.
    answer_signals = 3

    if calculated_signals == answer_signals:
        # The calculated number of signals matches the answer.
        return "Correct"
    else:
        # The calculation does not match the answer.
        reason = (
            f"Incorrect. The final product of the reaction sequence is 3,4-diethylhex-3-ene. "
            f"However, my code calculated that this molecule has {calculated_signals} 13C-NMR signals due to its symmetry. "
            f"The provided answer states there are {answer_signals} signals. The calculated value does not match the answer."
        )
        return reason

# To run the check, you would execute the following line:
# print(check_correctness())