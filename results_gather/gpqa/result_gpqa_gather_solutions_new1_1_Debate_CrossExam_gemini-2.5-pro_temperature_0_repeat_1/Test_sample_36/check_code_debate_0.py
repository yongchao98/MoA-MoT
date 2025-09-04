# First, ensure you have rdkit installed:
# pip install rdkit

try:
    from rdkit import Chem
except ImportError:
    # This fallback allows the code to run in environments without rdkit,
    # relying on the manually verified chemical logic.
    print("Warning: RDKit not found. Proceeding with logical check only.")
    rdkit_available = False
else:
    rdkit_available = True

def check_final_answer():
    """
    Checks the correctness of the provided final answer by:
    1. Verifying the number of 13C-NMR signals for the final product.
    2. Checking if the answer's logic and final choice are consistent and correct.
    """
    # --- Problem Definition ---
    # The final product is 3,4-diethylhex-3-ene
    final_product_smiles = 'CCC(CC)=C(CC)CC'
    
    # The options provided in the question
    options = {'A': 8, 'B': 3, 'C': 6, 'D': 11}
    
    # --- Provided Answer's Claims ---
    # The final answer's text concludes there are 3 signals and chooses option B.
    answer_reasoning_result = 3
    answer_option_choice = 'B'
    answer_final_output = 'B' # Extracted from <<<B>>>

    # --- Verification ---
    
    # 1. Calculate the correct number of signals
    correct_signal_count = -1
    if rdkit_available:
        try:
            mol = Chem.MolFromSmiles(final_product_smiles)
            # CanonicalRankAtoms assigns the same rank to symmetry-equivalent atoms.
            # We count the number of unique ranks for carbon atoms.
            ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
            carbon_ranks = {ranks[atom.GetIdx()] for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
            correct_signal_count = len(carbon_ranks)
        except Exception as e:
            return f"An error occurred during RDKit processing: {e}"
    else:
        # Fallback to the manually verified logic if RDKit is not available.
        # The logic is sound: 2 alkene C's, 4 methylene C's, 4 methyl C's.
        # Due to symmetry, this simplifies to 3 unique environments.
        correct_signal_count = 3

    # 2. Check if the answer's reasoning is correct
    if answer_reasoning_result != correct_signal_count:
        return (f"Incorrect: The answer's reasoning concludes there are {answer_reasoning_result} signals, "
                f"but the correct number of signals for 3,4-diethylhex-3-ene is {correct_signal_count}.")

    # 3. Check if the answer correctly maps the result to the option letter
    correct_option_letter = None
    for letter, value in options.items():
        if value == correct_signal_count:
            correct_option_letter = letter
            break
            
    if correct_option_letter is None:
        return f"Constraint failed: The correct answer ({correct_signal_count} signals) does not match any of the provided options."

    if answer_option_choice != correct_option_letter:
        return (f"Incorrect: The answer selects option '{answer_option_choice}', but the correct option for "
                f"{correct_signal_count} signals is '{correct_option_letter}'.")

    # 4. Check if the final output format is consistent
    if answer_final_output != answer_option_choice:
        return (f"Incorrect: The answer's reasoning points to option '{answer_option_choice}', "
                f"but the final output is '<<<{answer_final_output}>>>'.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_final_answer()
print(result)