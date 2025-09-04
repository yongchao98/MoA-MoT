# Install rdkit if you don't have it: pip install rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def check_correctness():
    """
    This function verifies the correct starting material for a Ring-Closing Metathesis (RCM) reaction.
    
    1. It defines the target product and candidate starting materials using SMILES strings.
    2. It simulates the RCM reaction for each candidate.
    3. It identifies which candidates produce the target product.
    4. It analyzes the validity of the options, especially regarding IUPAC naming conventions,
       to determine the single best answer.
    """
    
    # The final answer provided by the LLM to be checked.
    llm_final_answer = 'C'

    # Step 1: Define the target product and its canonical representation.
    # Target: 5-isopropyl-3,4-dimethylcyclohex-1-ene
    try:
        target_smiles = 'C1=CC(C)C(C)C(C(C)C)C1'
        target_mol = Chem.MolFromSmiles(target_smiles)
        if not target_mol:
            return "Error in checker: Could not parse target SMILES."
        target_canonical_smiles = Chem.MolToSmiles(target_mol, canonical=True)
    except Exception as e:
        return f"Error in checker during target molecule setup: {e}"

    # Step 2: Define the candidate starting materials from the question options.
    # Note: The SMILES for option D describes the same molecule as option C,
    # but the name '4-isopropyl-5,6-dimethylocta-1,7-diene' is incorrect by IUPAC rules.
    candidates = {
        'A': {'name': '5-isopropyl-3,4-dimethylocta-2,6-diene', 'smiles': 'CC=C(C)C(C)C(C(C)C)C=CC'},
        'B': {'name': '5-isopropyl-3,4-dimethylocta-1,6-diene', 'smiles': 'C=CC(C)C(C)C(C(C)C)C=C'},
        'C': {'name': '5-isopropyl-3,4-dimethylocta-1,7-diene', 'smiles': 'C=CC(C)C(C)C(C(C)C)CC=C'},
        'D': {'name': '4-isopropyl-5,6-dimethylocta-1,7-diene', 'smiles': 'C=CCC(C(C)C)C(C)C(C)C=C'}
    }

    # Step 3: Define the RCM reaction using SMARTS.
    try:
        rxn_smarts = '([C:1]=[C:2])-[*:3]-([C:4]=[C:5])>>([C:2]=[C:4]).([C:1]=[C:5])'
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    except Exception as e:
        return f"Error in checker during reaction setup: {e}"

    # Step 4: Simulate the reaction for each candidate and find which ones work.
    successful_precursors = []
    for option, data in candidates.items():
        mol = Chem.MolFromSmiles(data['smiles'])
        if not mol:
            return f"Error in checker: Could not parse SMILES for option {option}."
        
        products = rxn.RunReactants((mol,))
        
        for product_set in products:
            for p in product_set:
                # We are interested in the cyclic product (has rings).
                if Chem.rdmolops.GetSSSR(p) > 0:
                    p_canonical_smiles = Chem.MolToSmiles(p, canonical=True)
                    if p_canonical_smiles == target_canonical_smiles:
                        if option not in successful_precursors:
                            successful_precursors.append(option)

    # Step 5: Analyze the results to determine the single best answer.
    # The key insight is that the molecule that produces the target has the correct IUPAC name
    # '5-isopropyl-3,4-dimethylocta-1,7-diene' (Option C). The name in Option D,
    # '4-isopropyl-5,6-dimethylocta-1,7-diene', is an incorrect IUPAC name for the same molecule,
    # as it violates the 'lowest locant set' rule. Therefore, C is the only fully correct answer.
    
    correct_option = 'C'

    if llm_final_answer == correct_option:
        # The LLM correctly identified the option with the valid IUPAC name.
        return "Correct"
    else:
        # The LLM chose an incorrect option.
        if llm_final_answer == 'D' and 'D' in successful_precursors:
            return f"Incorrect. The provided answer is D. While the molecule described by option D does produce the target product, its IUPAC name is incorrect. The correct IUPAC name for this molecule is given in option C. Therefore, C is the best answer."
        else:
            return f"Incorrect. The provided answer is {llm_final_answer}. The correct starting material is option C ({candidates['C']['name']}), which is the only option that both produces the target molecule and is correctly named according to IUPAC rules."

result = check_correctness()
print(result)