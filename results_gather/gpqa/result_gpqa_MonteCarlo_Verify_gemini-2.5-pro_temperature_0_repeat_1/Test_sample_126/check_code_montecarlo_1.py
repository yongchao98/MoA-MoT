try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import rdBase
    # Suppress RDKit's verbose logging for a cleaner output
    rdBase.DisableLog('rdApp.warning')
except ImportError:
    # Handle cases where the required library is not installed.
    print("RDKit library not found. Please install it (e.g., 'pip install rdkit-pypi') to run this chemical validation.")
    exit()

def check_cope_rearrangement():
    '''
    This function checks the product of a Cope rearrangement on 5-butylnona-2,6-diene.
    1. It defines the starting molecule using its SMILES string.
    2. It defines the Cope rearrangement reaction using a reaction SMARTS pattern.
    3. It runs the reaction to generate the product(s).
    4. It converts the product molecule to a canonical SMILES string for unambiguous comparison.
    5. It compares the product's SMILES with the SMILES of the provided options.
    '''
    
    # --- Step 1: Define Molecules ---
    
    # Starting molecule: 5-butylnona-2,6-diene
    # SMILES representation: CC=CCC(CCCC)C=CCC
    start_smiles = "CC=CCC(CCCC)C=CCC"
    start_mol = Chem.MolFromSmiles(start_smiles)
    
    if start_mol is None:
        return "Error: Could not parse the starting molecule's SMILES string."

    # The LLM's answer to be checked
    llm_answer_key = "A"

    # Options provided in the question, represented by their canonical SMILES.
    # This avoids ambiguity from different valid SMILES strings for the same molecule.
    options = {
        "A": "CCCCCC=C(CC)C(C)C=C",  # Canonical SMILES for 4-ethyl-3-methyldeca-1,5-diene
        "B": "CCCC=C(CC)C(C)C=CC",  # Canonical SMILES for 5-ethyl-4-methyldeca-2,6-diene
        "C": "CCCCC=C(CC)CCC=CC",   # Canonical SMILES for 5-ethylundeca-2,6-diene
        "D": "CCCC=C(CC)C(C)C=CC"   # Same as B
    }

    # --- Step 2: Define the Reaction ---

    # A Cope rearrangement is a [3,3]-sigmatropic rearrangement.
    # The reaction SMARTS describes the transformation of a hexa-1,5-diene system.
    # Bonds broken: C3-C4. Bonds made: C1-C6. Double bonds shift.
    # SMARTS: [C:1]=[C:2]-[C:3]-[C:4]-[C:5]=[C:6] >> [C:2]=[C:3]-[C:4]=[C:5]-[C:6]-[C:1]
    # Note: Atom maps in RDKit start from 1.
    rxn_smarts = '[C:1]=[C:2][C:3][C:4][C:5]=[C:6]>>[C:2]=[C:3][C:4]=[C:5][C:6][C:1]'
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    except Exception as e:
        return f"Error: Could not create reaction from SMARTS pattern. Details: {e}"

    # --- Step 3: Run the Reaction ---

    products = rxn.RunReactants((start_mol,))

    if not products:
        return "Failure: The reaction did not yield any products. The SMARTS pattern may not have matched the starting molecule."

    # The reaction returns a tuple of tuples of product molecules. We expect one product set.
    product_mol = products[0][0]
    
    try:
        # Sanitize the molecule to ensure correct valences and properties.
        Chem.SanitizeMol(product_mol)
    except Exception as e:
        return f"Error: Could not sanitize the product molecule. Details: {e}"

    # --- Step 4: Compare and Verify ---

    # Generate the canonical SMILES for the actual product (ignoring stereochemistry for this problem).
    actual_product_smiles = Chem.MolToSmiles(product_mol, isomericSmiles=False, canonical=True)
    
    # Get the canonical SMILES for the LLM's answer.
    llm_answer_smiles = options.get(llm_answer_key)

    if actual_product_smiles == llm_answer_smiles:
        return "Correct"
    else:
        # If the answer is wrong, find the correct one.
        correct_option_key = None
        for key, smiles in options.items():
            if smiles == actual_product_smiles:
                correct_option_key = key
                break
        
        reason = f"Incorrect. The provided answer was {llm_answer_key}.\n"
        reason += f"The starting molecule 5-butylnona-2,6-diene undergoes a Cope rearrangement.\n"
        reason += f"The calculated product is 4-ethyl-3-methyldeca-1,5-diene, which has the canonical SMILES: {actual_product_smiles}.\n"
        reason += f"The molecule for option {llm_answer_key} has the SMILES: {llm_answer_smiles}.\n"
        if correct_option_key:
            reason += f"The correct option is {correct_option_key}."
        else:
            reason += "None of the provided options match the correct product."
            
        return reason

result = check_cope_rearrangement()
print(result)