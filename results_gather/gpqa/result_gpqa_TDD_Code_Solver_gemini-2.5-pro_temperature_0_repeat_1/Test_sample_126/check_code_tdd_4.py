import sys

def run_check():
    """
    This function contains the Python code to check the correctness of the chemical reaction answer.
    It will be executed to print the verification result.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
    except ImportError:
        print("Skipped: The RDKit library is required for this chemical verification but is not installed. Please install it using 'pip install rdkit-pypi'.")
        return

    def get_molecule_info_from_smiles(name, smiles_string):
        """
        Helper function to create an RDKit molecule object from a SMILES string
        and calculate its molecular formula and canonical SMILES.
        """
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return None, None, None
        
        # Add explicit hydrogens to get the correct formula
        mol_with_h = Chem.AddHs(mol)
        formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol_with_h)
        
        # Generate a canonical SMILES string for consistent comparison.
        # We use isomericSmiles=False because we are checking for constitutional isomers,
        # and stereochemistry is not specified in the question.
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
        
        return mol, formula, canonical_smiles

    def check_cope_rearrangement():
        """
        Checks if the Cope rearrangement of 5-butylnona-2,6-diene yields
        4-ethyl-3-methyldeca-1,5-diene.
        """
        # --- Step 1: Define Reactant and Proposed Product ---
        # IUPAC names and their corresponding SMILES strings.
        # Using manually verified SMILES is more reliable than name-to-structure converters for complex molecules.
        reactant_name = "5-butylnona-2,6-diene"
        reactant_smiles = "CC=CCC(CCCC)C=CCC"

        proposed_answer_name = "4-ethyl-3-methyldeca-1,5-diene" # This is option B from the question
        proposed_answer_smiles = "C=CC(C)(C(CC)C=CCCCC)" # A more canonical representation of C=CC(C)C(CC)C=CCCCC

        # --- Step 2: Get Molecular Representations ---
        reactant_mol, reactant_formula, _ = get_molecule_info_from_smiles(reactant_name, reactant_smiles)
        product_mol_B, product_formula_B, product_smiles_B = get_molecule_info_from_smiles(proposed_answer_name, proposed_answer_smiles)

        if not reactant_mol or not product_mol_B:
            return "Error: Failed to parse molecule SMILES strings. Cannot verify."

        # --- Step 3: Check Basic Constraints ---
        # Constraint 1: Isomerism. A rearrangement reaction must produce an isomer.
        if reactant_formula != product_formula_B:
            return f"Incorrect. The product must be an isomer of the reactant. Reactant formula is {reactant_formula}, but the proposed product (B) has formula {product_formula_B}."

        # --- Step 4: Simulate the Reaction ---
        # The reaction is a [3,3]-sigmatropic Cope rearrangement.
        # We model this with a reaction SMARTS string.
        # The mapping is: [C:1]=[C:2]-[C:3]-[C:4]-[C:5]=[C:6]>>[C:2]=[C:3]-[C:4]=[C:5]-[C:6]-[C:1]
        # This represents the breaking of the C3-C4 bond, formation of the C1-C6 bond, and shift of the pi bonds.
        try:
            rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]-[C:3]-[C:4]-[C:5]=[C:6]>>[C:2]=[C:3]-[C:4]=[C:5]-[C:6]-[C:1]')
        except Exception as e:
            return f"Error: Could not define the reaction in RDKit: {e}"

        # Constraint 2: Reactant Suitability. Check if the reactant has the required 1,5-diene substructure.
        if not reactant_mol.HasSubstructMatch(rxn.GetReactantTemplate(0)):
            return f"Error: The reactant '{reactant_name}' does not contain the required 1,5-diene substructure for the defined Cope rearrangement. The premise of the question might be flawed or the reactant SMILES is incorrect."

        # Run the reaction on the reactant molecule
        products = rxn.RunReactants((reactant_mol,))

        if not products or not products[0]:
            return "Incorrect. The simulated Cope rearrangement did not yield any products for the given reactant. While the structure is suitable, the simulation failed."

        # The reaction yields a tuple of tuples of product molecules. We expect one product set.
        reaction_product_mol = products[0][0]
        
        # Sanitize the product molecule to ensure correct valences and structure
        try:
            Chem.SanitizeMol(reaction_product_mol)
        except Exception as e:
            return f"Error: Could not sanitize the reaction product molecule: {e}"

        # Get the canonical, non-isomeric SMILES of the simulated product
        simulated_product_smiles = Chem.MolToSmiles(reaction_product_mol, isomericSmiles=False, canonical=True)

        # --- Step 5: Compare Simulated Product with Proposed Answer ---
        if simulated_product_smiles == product_smiles_B:
            return "Correct"
        else:
            return f"Incorrect. The simulated product of the Cope rearrangement does not match the proposed answer (B).\\nReactant: {reactant_name} ({reactant_smiles})\\nSimulated Product SMILES: {simulated_product_smiles}\\nProposed Answer B SMILES: {product_smiles_B}"

    # --- Execute the check and print the result ---
    result = check_cope_rearrangement()
    print(result)

# To run this check, you need to have the RDKit library installed.
# You can install it via pip:
# pip install rdkit-pypi
run_check()