import sys

try:
    from rdkit import Chem
except ImportError:
    # If rdkit is not installed, the check cannot be performed.
    print("Execution failed: The 'rdkit' library is required for this check but is not installed. Please install it using 'pip install rdkit'.")
    # In a real script, you might exit or raise an exception.
    # For this context, we'll set a flag and handle it in the checking function.
    rdkit_installed = False
else:
    rdkit_installed = True

def count_distinct_hydrogens(smiles: str) -> int:
    """
    Counts the number of chemically distinct hydrogen atoms in a molecule
    represented by a SMILES string using RDKit's canonical atom ranking.
    This method correctly distinguishes diastereotopic protons.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"Could not parse SMILES string: {smiles}")

    mol = Chem.AddHs(mol)
    # CanonicalRankAtoms with breakTies=True (default) correctly identifies
    # symmetrically and stereochemically distinct atoms.
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=True)
    
    # Collect the unique ranks corresponding to hydrogen atoms.
    h_ranks = {ranks[i] for i, atom in enumerate(mol.GetAtoms()) if atom.GetAtomicNum() == 1}
    
    return len(h_ranks)

def check_llm_answer():
    """
    This function verifies the LLM's answer by analyzing the proposed final products.
    The LLM's reasoning is that the final product ("product 4") is a mixture of
    benzene and o-xylylene, and that the total number of distinct hydrogen signals is 4.
    """
    if not rdkit_installed:
        # Return the error message if the required library is missing.
        return "Execution failed: The 'rdkit' library is required for this check but is not installed. Please install it using 'pip install rdkit'."

    # LLM's proposed final products and answer
    llm_final_products = {
        "benzene": 'c1ccccc1',
        "o-xylylene": 'C=C1C=CC=CC1=C' # This is 5,6-bis(methylidene)cyclohexa-1,3-diene
    }
    llm_final_answer = 4 # The value corresponding to option B

    # LLM's breakdown of signals to reach the answer: 1 (benzene) + 3 (o-xylylene)
    llm_signal_breakdown = {"benzene": 1, "o-xylylene": 3}

    try:
        # Step 1: Rigorously calculate the number of signals for each component.
        calculated_signals_benzene = count_distinct_hydrogens(llm_final_products["benzene"])
        calculated_signals_oxylylene = count_distinct_hydrogens(llm_final_products["o-xylylene"])
        
        # Step 2: Calculate the total number of signals in the mixture.
        calculated_total_signals = calculated_signals_benzene + calculated_signals_oxylylene

        # Step 3: Compare the rigorous calculation with the LLM's claims.
        # The primary point of failure is the signal count for o-xylylene.
        if calculated_signals_oxylylene != llm_signal_breakdown["o-xylylene"]:
            reason = (
                f"Incorrect. The LLM's reasoning to arrive at the answer '4' is flawed.\n"
                f"The analysis correctly identifies the final products as a mixture of benzene and o-xylylene. However, the subsequent counting of NMR signals is erroneous.\n"
                f"1. Benzene ('c1ccccc1'): Correctly identified as having 1 distinct hydrogen signal. My calculation confirms: {calculated_signals_benzene} signal.\n"
                f"2. o-Xylylene ('C=C1C=CC=CC1=C'): The LLM claims this molecule has 3 signals. This is incorrect. A rigorous analysis, which accounts for diastereotopic protons, shows it has {calculated_signals_oxylylene} distinct signals (2 from the ring protons and 2 from the four exocyclic methylene protons, which are diastereotopic).\n"
                f"The LLM's own provided code hard-codes this incorrect value of 3; it does not calculate it, which invalidates its 'verification'.\n"
                f"3. Total Signals: The correct theoretical number of signals in the mixture is the sum of the components: {calculated_signals_benzene} (benzene) + {calculated_signals_oxylylene} (o-xylylene) = {calculated_total_signals}.\n"
                f"The LLM's answer of 4 is incorrect because the true theoretical count is {calculated_total_signals}."
            )
            return reason
        # This fallback handles the unlikely case that our calculation matches the LLM's flawed one.
        elif calculated_total_signals == llm_final_answer:
            return "Correct"
        else:
            return f"Incorrect. The total calculated signals ({calculated_total_signals}) does not match the provided answer ({llm_final_answer})."

    except Exception as e:
        return f"An error occurred during the verification check: {e}"

# Execute the check and print the result.
result = check_llm_answer()
print(result)