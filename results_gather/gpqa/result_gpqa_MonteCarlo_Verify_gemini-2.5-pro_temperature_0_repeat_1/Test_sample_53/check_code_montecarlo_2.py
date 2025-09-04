import sys

try:
    # RDKit is a standard cheminformatics library for processing chemical structures.
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    # The check cannot run without RDKit.
    print("Error: RDKit library not found. Please install it (`pip install rdkit-pypi`) to run this verification code.")
    sys.exit(1)

def check_rcm_correctness():
    """
    Checks if the starting material from the given answer correctly forms the target 
    product via Ring-Closing Metathesis (RCM).
    """
    # --- Problem Definition ---
    # Target Product: 5-isopropyl-3,4-dimethylcyclohex-1-ene
    # Reaction: Ring-Closing Metathesis (RCM)
    # Provided Answer to Check: D

    # --- Data Setup ---
    # The canonical SMILES for the target product is used as the ground truth.
    # This was obtained from chemical databases (e.g., PubChem CID 535319).
    target_canonical_smiles = "CC1C(C)C(C(C)C)CC=C1"

    # The option selected by the LLM to be checked.
    selected_option = "D"
    
    # Data for the selected starting material. The SMILES string is derived from its IUPAC name.
    reactant_data = {
        "name": "5-isopropyl-3,4-dimethylocta-1,7-diene",
        "smiles": "C=CC(C)C(C)C(C(C)C)CC=C"
    }

    # --- Chemical Simulation ---
    
    reactant_smiles = reactant_data["smiles"]
    reactant_name = reactant_data["name"]

    # Create an RDKit molecule object from the SMILES string.
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    if not reactant_mol:
        return f"Error: Could not parse the SMILES string for option {selected_option}: '{reactant_smiles}'"

    # Define the RCM reaction using a reaction SMARTS string. This is a general
    # definition for metathesis that RDKit applies intramolecularly for dienes.
    # It finds two C=C bonds and swaps their ends, forming the ring and a byproduct.
    try:
        rxn = AllChem.ReactionFromSmarts('([C:1]=[C:2]).([C:3]=[C:4])>>[C:2]=[C:3].[C:1]=[C:4]')
        products = rxn.RunReactants((reactant_mol,))
    except Exception as e:
        return f"An error occurred during the RDKit reaction simulation: {e}"

    # --- Analysis of Results ---

    if not products:
        return (f"Incorrect. The starting material for option {selected_option} ({reactant_name}) "
                f"did not yield any products under the simulated RCM reaction, "
                f"suggesting it is not a suitable substrate.")

    # The result is a tuple of product sets. We expect one set for one reactant.
    product_set = products[0]
    
    # Identify the major product (the cyclic compound) by its size (atom count).
    # The byproduct (ethene, C2H4) will have fewer atoms.
    main_product_mol = None
    if len(product_set) > 0:
        sorted_products = sorted(list(product_set), key=lambda m: m.GetNumHeavyAtoms(), reverse=True)
        main_product_mol = sorted_products[0]

    if not main_product_mol:
        return (f"Incorrect. The RCM reaction for option {selected_option} did not produce an identifiable main product.")

    # Sanitize the resulting molecule to ensure correct valencies and structure.
    try:
        Chem.SanitizeMol(main_product_mol)
    except Exception:
        return (f"Incorrect. The RCM reaction for option {selected_option} produced a chemically invalid structure.")

    # Generate the canonical SMILES for the simulated product for a standard comparison.
    simulated_product_canonical_smiles = Chem.MolToSmiles(main_product_mol, canonical=True)

    # --- Final Verdict ---
    # Compare the SMILES of the simulated product with the target product.
    if simulated_product_canonical_smiles == target_canonical_smiles:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {selected_option} ({reactant_name}). "
                f"Ring-closing metathesis of this material produces a product with canonical SMILES '{simulated_product_canonical_smiles}'. "
                f"The target product, 5-isopropyl-3,4-dimethylcyclohex-1-ene, has canonical SMILES '{target_canonical_smiles}'. "
                f"Since the SMILES strings do not match, the answer is incorrect.")

# Execute the check and print the result.
result = check_rcm_correctness()
print(result)