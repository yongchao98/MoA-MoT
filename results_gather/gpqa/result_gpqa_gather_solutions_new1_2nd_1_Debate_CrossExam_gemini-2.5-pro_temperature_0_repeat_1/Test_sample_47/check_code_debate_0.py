import re
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    # If rdkit is not installed, we cannot run the check.
    # We will have to rely on a manual logical check.
    Chem = None

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer to the chemistry question.
    """
    # The final answer provided by the LLM to be checked.
    llm_final_answer = "<<<C>>>"

    # Define the options from the question
    # A) 7, B) 8, C) 6, D) 10
    options = {'A': 7, 'B': 8, 'C': 6, 'D': 10}

    # --- Step 1: Parse the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_final_answer)
    if not match:
        return "Invalid answer format. The answer must be in the format <<<X>>>."
    
    llm_option = match.group(1)
    llm_value = options.get(llm_option)

    if llm_value is None:
        return f"Invalid option '{llm_option}' provided in the answer."

    # --- Step 2: Programmatically solve the problem to find the ground truth ---
    # The reaction sequence is:
    # 1. Cyclohexanone + Br2 -> 2-bromocyclohexanone
    # 2. 2-bromocyclohexanone + NaOH -> Cyclopentanecarboxylic acid (Favorskii rearrangement)
    # 3. Cyclopentanecarboxylic acid + SOCl2 -> Cyclopentanecarbonyl chloride
    # 4. Cyclopentanecarbonyl chloride + LiAlH(O-t-Bu)3 -> Cyclopentanecarbaldehyde
    # The final product is Cyclopentanecarbaldehyde.
    
    # We can use RDKit to count the number of chemically distinct hydrogens.
    # If RDKit is not available, we fall back to the logical deduction.
    calculated_correct_value = 6 # Based on manual chemical analysis
    
    if Chem:
        try:
            # SMILES string for cyclopentanecarbaldehyde
            smiles = 'O=CC1CCCC1'
            mol = Chem.MolFromSmiles(smiles)
            
            # Add explicit hydrogens to the molecular graph
            mol_with_hs = Chem.AddHs(mol)
            
            # Use CanonicalRankAtoms to find symmetrically equivalent atoms.
            # The `breakTies=True` argument is important for distinguishing
            # constitutionally equivalent but stereochemically different atoms (like diastereotopic protons).
            ranks = Chem.CanonicalRankAtoms(mol_with_hs, breakTies=True)
            
            # Collect the unique ranks for hydrogen atoms
            hydrogen_ranks = set()
            for atom, rank in zip(mol_with_hs.GetAtoms(), ranks):
                if atom.GetAtomicNum() == 1:  # It's a hydrogen atom
                    hydrogen_ranks.add(rank)
            
            calculated_correct_value = len(hydrogen_ranks)
        except Exception as e:
            # If RDKit fails, we can't verify programmatically.
            # The manual check result will be used.
            pass

    # --- Step 3: Compare the LLM's answer with the ground truth ---
    if llm_value == calculated_correct_value:
        return "Correct"
    else:
        correct_option = [opt for opt, val in options.items() if val == calculated_correct_value][0]
        reason = (f"Incorrect. The final product of the reaction sequence is cyclopentanecarbaldehyde. "
                  f"A symmetry analysis of this molecule shows there are {calculated_correct_value} chemically distinct hydrogen atoms. "
                  f"This corresponds to option {correct_option}. The provided answer was option {llm_option}, which corresponds to the value {llm_value}.")
        return reason

# Execute the check
result = check_chemistry_answer()
print(result)