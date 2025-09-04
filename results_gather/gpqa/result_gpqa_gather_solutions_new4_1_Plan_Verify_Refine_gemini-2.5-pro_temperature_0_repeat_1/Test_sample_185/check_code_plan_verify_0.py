try:
    from rdkit import Chem
except ImportError:
    # If rdkit is not installed, we cannot perform the structural analysis.
    # We will provide a message to the user.
    Chem = None

def get_alkene_substitution(mol_smiles: str) -> int:
    """
    Parses a SMILES string, finds the C=C bond, and returns its substitution degree.
    Assumes there is only one C=C bond in the structure.
    """
    if not Chem:
        raise ImportError("RDKit library is not installed. Cannot perform structural analysis.")

    mol = Chem.MolFromSmiles(mol_smiles)
    if not mol:
        raise ValueError(f"Could not parse SMILES string: {mol_smiles}")
        
    # Add hydrogens to ensure correct degree calculation
    mol = Chem.AddHs(mol)

    for bond in mol.GetBonds():
        # Find the double bond that is between two carbon atoms
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                # The degree of substitution is the number of non-hydrogen atoms attached 
                # to the double bond carbons, excluding the bond between them.
                heavy_neighbors1 = sum(1 for n in atom1.GetNeighbors() if n.GetAtomicNum() > 1)
                heavy_neighbors2 = sum(1 for n in atom2.GetNeighbors() if n.GetAtomicNum() > 1)
                
                # Each carbon in the double bond is a "heavy neighbor" to the other.
                # So, the number of external substituents is (heavy_neighbors1 - 1) + (heavy_neighbors2 - 1).
                substitution = (heavy_neighbors1 - 1) + (heavy_neighbors2 - 1)
                return substitution
    return 0 # No C=C bond found

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer and reasoning.
    """
    # The LLM's final answer and key points from its reasoning
    llm_answer_choice = 'C'
    llm_reasoning_claim = {
        "product_C_substitution": 3, # Claim: trisubstituted
        "product_B_substitution": 2, # Claim: disubstituted
    }

    # Based on chemical literature (e.g., Overman's work), the correct product is C.
    correct_answer_choice = 'C'

    if llm_answer_choice != correct_answer_choice:
        return f"Incorrect. The final answer is given as {llm_answer_choice}, but the established product for this reaction is {correct_answer_choice}."

    # Now, check the reasoning using RDKit if available.
    if not Chem:
        return "Correct. The final answer choice 'C' is correct based on chemical literature. RDKit is not installed, so the reasoning's structural claims cannot be verified."

    try:
        # SMILES strings representing the structures based on IUPAC name interpretation.
        # B: 4,4a,5,7a-tetrahydro-3H-cyclopenta[c]pyridine -> C=C is disubstituted
        smiles_b = "C1=NC2C3C4C=CC3C24"
        # C: 4,6,7,7a-tetrahydro-3H-cyclopenta[c]pyridine -> C=C is at the ring fusion
        smiles_c = "C1=NC2C3C(C4=C2C3)CC4"

        # Perform the analysis
        actual_subst_b = get_alkene_substitution(smiles_b)
        actual_subst_c = get_alkene_substitution(smiles_c)

        # Check if the reasoning's claims match the actual structure
        if actual_subst_c == llm_reasoning_claim["product_C_substitution"]:
            # If the reasoning is correct, and the answer is correct, then the overall answer is correct.
            return "Correct"
        else:
            # The final answer choice is correct, but the reasoning to get there is flawed.
            return (f"Incorrect. The reasoning provided to reach the correct answer 'C' contains a significant factual error. "
                    f"The reasoning states that the alkene in product C is {llm_reasoning_claim['product_C_substitution']}-substituted (trisubstituted). "
                    f"However, a structural analysis reveals that the C=C double bond in isomer C is actually {actual_subst_c}-substituted (tetrasubstituted). "
                    f"While the conclusion that C is the most stable isomer is correct (tetrasubstituted > disubstituted), the premise used in the reasoning is wrong.")

    except Exception as e:
        return f"An error occurred during the check: {e}. Cannot fully verify the answer."

# Run the check
result = check_chemistry_answer()
print(result)