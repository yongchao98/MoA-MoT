from rdkit import Chem
from rdkit.Chem import AllChem

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer by:
    1. Defining the final product based on the reaction sequence.
    2. Using RDKit to analyze the product's structure and symmetry.
    3. Calculating the number of unique carbon environments (13C-NMR signals).
    4. Comparing the calculated result with the provided answer's reasoning and final choice.
    """
    try:
        # --- Step 1: Define the problem and the provided answer ---
        
        # Options from the question
        options = {'A': 6, 'B': 3, 'C': 8, 'D': 11}
        
        # The final answer from the LLM being checked
        llm_reasoning_value = 3
        llm_final_choice = 'B'

        # --- Step 2: Determine the final product structure ---
        
        # The reaction sequence is:
        # 1. Propionaldehyde -> 2-ethyl-1,3-dithiolane
        # 2. Deprotonation -> Lithiated carbanion
        # 3. Alkylation with bromoethane -> 2,2-diethyl-1,3-dithiolane
        # 4. Deprotection -> 3-pentanone (diethyl ketone)
        # 5. Wittig with ylide from 3-bromopentane -> 3,4-diethylhex-3-ene
        # The SMILES string for 3,4-diethylhex-3-ene is CCC(=C(CC)CC)CC
        final_product_smiles = "CCC(=C(CC)CC)CC"

        # --- Step 3: Analyze the final product using RDKit ---
        
        mol = Chem.MolFromSmiles(final_product_smiles)
        if mol is None:
            return "Error in checker: Could not create molecule from SMILES string."

        # A check to confirm that E/Z isomerism is not possible for this molecule.
        # EnumerateStereoisomers will return only one molecule if no stereoisomers exist.
        isomers = list(Chem.EnumerateStereoisomers.EnumerateStereoisomers(mol))
        if len(isomers) > 1:
            return f"Incorrect. The analysis missed that the final product has {len(isomers)} stereoisomers, which would complicate the NMR spectrum. The provided answer's reasoning about symmetry is flawed."

        # Add hydrogens to ensure the molecular graph is complete for symmetry perception
        mol = Chem.AddHs(mol)

        # Use CanonicalRankAtoms to find symmetrically equivalent atoms.
        # Atoms with the same rank are considered equivalent.
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
        
        # Filter for carbon atoms (atomic number 6) and find their unique ranks
        carbon_ranks = set()
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                carbon_ranks.add(ranks[atom.GetIdx()])
        
        calculated_signals = len(carbon_ranks)

        # --- Step 4: Verify the provided answer ---

        # Check 1: Does the calculated number of signals match the reasoning?
        if calculated_signals != llm_reasoning_value:
            return f"Incorrect. The reasoning states there are {llm_reasoning_value} signals, but a computational analysis of 3,4-diethylhex-3-ene shows {calculated_signals} signals."

        # Check 2: Does the chosen option letter match the value from the reasoning?
        if options.get(llm_final_choice) != llm_reasoning_value:
            return f"Incorrect. The reasoning correctly identifies {llm_reasoning_value} signals, but the final choice '{llm_final_choice}' corresponds to the value {options.get(llm_final_choice)} in the options list, which is a mismatch."

        # If all checks pass, the answer is correct.
        return "Correct"

    except ImportError:
        return "Error: RDKit library not found. Please install it using 'pip install rdkit'."
    except Exception as e:
        return f"An unexpected error occurred: {e}"

# Execute the check
result = check_answer_correctness()
print(result)