import sys
# It's good practice to check for the library and provide a helpful message if it's not there.
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    # If rdkit is not installed, we cannot perform the check.
    # We will print an error message to stderr and exit, so the calling environment knows something went wrong.
    sys.stderr.write("Error: RDKit library not found. Please install it (`pip install rdkit`) to run the checker.\n")
    # Returning a specific string to indicate the failure.
    print("Failure: RDKit not installed.")
    # sys.exit(1) # In a real script, you might exit. Here we'll just return a message.

def check_correctness():
    '''
    Checks if the starting material from option A correctly synthesizes 
    5-isopropyl-3,4-dimethylcyclohex-1-ene via Ring-Closing Metathesis (RCM).
    '''
    # Define the target product by its structure (SMILES) and get its canonical form for unique comparison.
    target_name = "5-isopropyl-3,4-dimethylcyclohex-1-ene"
    # The SMILES for the target product is derived from its IUPAC name: a cyclohexene ring with
    # methyl groups at C3 and C4, and an isopropyl group at C5.
    target_smiles = "CC(C)C1CC(C)C(C)C=C1" 
    target_mol = Chem.MolFromSmiles(target_smiles)
    # A canonical SMILES is a standardized representation, essential for accurate string comparison of molecules.
    target_canonical_smiles = Chem.MolToSmiles(target_mol, canonical=True)

    # Define the proposed starting material from the provided answer (A).
    proposed_answer_key = "A"
    options = {
        "A": "4-isopropyl-5,6-dimethylocta-1,7-diene",
        "B": "5-isopropyl-3,4-dimethylocta-2,6-diene",
        "C": "5-isopropyl-3,4-dimethylocta-1,7-diene",
        "D": "5-isopropyl-3,4-dimethylocta-1,6-diene"
    }
    smiles_dict = {
        "A": "C=CCC(C(C)C)C(C)C(C)C=C",
        "B": "CC=CC(C)C(C(C)C)C=CC",
        "C": "C=CC(C)C(C)C(C(C)C)CC=C",
        "D": "C=CC(C)C(C)C(C(C)C)C=CC"
    }
    proposed_smiles = smiles_dict[proposed_answer_key]
    reactant_mol = Chem.MolFromSmiles(proposed_smiles)

    # --- Constraint Checking ---
    # Constraint 1: The reaction is RCM to form a 6-membered ring. This requires a 1,7-diene.
    # Options B (2,6-diene) and D (1,6-diene) are incorrect based on this constraint.
    # The proposed answer A is a 1,7-diene, so it passes this constraint.

    # --- Reaction Simulation ---
    # Define the RCM reaction using a reaction SMARTS string.
    # This SMARTS will find any 1,7-diene chain and cyclize it, forming a cyclohexene and ethene.
    rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]-[*:3]-[*:4]-[*:5]-[*:6]-[C:7]=[C:8]>>[C:2]1-[*:3]-[*:4]-[*:5]-[*:6]-[C:7]=1.[C:1]=[C:8]')
    
    # Run the reaction on the proposed starting material.
    products = rxn.RunReactants((reactant_mol,))

    # Check if the reaction produced anything. If not, the reactant didn't match the reaction pattern.
    if not products:
        return f"Incorrect. The proposed starting material, {options[proposed_answer_key]}, does not match the required 1,7-diene structure for this RCM reaction."

    # The reaction yields a tuple of product sets. We are interested in the first product of the first set (the cycloalkene).
    product_mol = products[0][0]
    Chem.SanitizeMol(product_mol)
    product_canonical_smiles = Chem.MolToSmiles(product_mol, canonical=True)

    # --- Verification ---
    # Compare the canonical SMILES of the simulated product with the target's canonical SMILES.
    if product_canonical_smiles == target_canonical_smiles:
        return "Correct"
    else:
        # If they don't match, the answer is incorrect. We can determine the name of the actual product for a detailed explanation.
        # The product from option A is not the target. Let's find the name of the isomer that was actually formed.
        # The product from option C, 5-isopropyl-3,4-dimethylocta-1,7-diene, actually yields 4-isopropyl-5,6-dimethylcyclohex-1-ene.
        # The product from option A, 4-isopropyl-5,6-dimethylocta-1,7-diene, yields the target.
        # Therefore, this 'else' block should not be reached if the provided answer is A.
        return (f"Incorrect. The RCM product of the proposed starting material ({options[proposed_answer_key]}) "
                f"has the canonical SMILES '{product_canonical_smiles}', which does not match the target's "
                f"canonical SMILES '{target_canonical_smiles}'. The actual product formed is an isomer of the target.")

# Execute the check function and print the result.
# This structure ensures that the check is only run when the script is executed directly.
if 'rdkit' in sys.modules:
    result = check_correctness()
    print(result)