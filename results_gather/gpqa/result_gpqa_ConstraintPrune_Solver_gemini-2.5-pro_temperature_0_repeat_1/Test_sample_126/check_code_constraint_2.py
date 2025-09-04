# The RDKit library is required to run this code.
# You can install it with: pip install rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("RDKit not found. Please install it using 'pip install rdkit' to run this check.")
    # Define a dummy function to avoid errors if rdkit is not installed.
    def check_reaction_correctness():
        return "Skipped: RDKit library is not installed. This check cannot be performed."
else:
    def check_reaction_correctness():
        """
        Checks the correctness of the LLM's answer for the Cope rearrangement of 5-butylnona-2,6-diene.
        It simulates the reaction using RDKit and compares the product to the provided options.
        """
        try:
            # Step 1: Define the reactant and candidate molecules using SMILES.
            # These are taken from the LLM's answer to verify its reasoning.
            structures = {
                "Reactant": {"name": "5-butylnona-2,6-diene", "smiles": "CCC=CC(CCCC)C=CCC"},
                "A": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCC=C(CC)C(C)C=C"},
                "B": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=C(C)C(CC)C=CC"},
                "C": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCCC=C(CC)C=CCC"},
                "D": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=C(C)C(CC)C=CC"}
            }
            llm_answer_key = "A"

            # Generate canonical SMILES for all structures for unambiguous comparison.
            canonical_smiles = {}
            for key, props in structures.items():
                mol = Chem.MolFromSmiles(props["smiles"])
                if mol is None:
                    return f"Error: The SMILES string for '{props['name']}' ('{props['smiles']}') is invalid."
                # Using isomericSmiles=True is good practice for stereochemically-aware comparison.
                canonical_smiles[key] = Chem.MolToSmiles(mol, isomericSmiles=True)

            reactant_mol = Chem.MolFromSmiles(structures["Reactant"]["smiles"])

            # Step 2: Define the Cope rearrangement reaction using SMARTS.
            # The SMARTS pattern correctly maps the [3,3]-sigmatropic shift.
            rxn_smarts = '[C:1]=[C:2]-[C:3]-[C:4]-[C:5]=[C:6]>>[C:2]=[C:3]-[C:4]=[C:5]-[C:6]-[C:1]'
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            
            # Step 3: Run the reaction simulation.
            products = rxn.RunReactants((reactant_mol,))

            if not products or not products[0]:
                return "Incorrect: The reaction simulation failed to produce any products from the reactant."

            # Step 4: Get the canonical SMILES of the simulated product.
            product_mol = products[0][0]
            Chem.SanitizeMol(product_mol)
            product_canonical_smiles = Chem.MolToSmiles(product_mol, isomericSmiles=True)

            # Step 5: Compare the simulated product with the LLM's chosen answer.
            llm_answer_smiles = canonical_smiles[llm_answer_key]

            if product_canonical_smiles == llm_answer_smiles:
                return "Correct"
            else:
                # If the answer is wrong, identify the correct option for a more informative message.
                correct_key = None
                for key, smiles in canonical_smiles.items():
                    if key != "Reactant" and smiles == product_canonical_smiles:
                        correct_key = key
                        break
                
                if correct_key:
                    return f"Incorrect: The LLM's answer is {llm_answer_key}, but the simulation shows the correct product is {correct_key} ({structures[correct_key]['name']})."
                else:
                    return f"Incorrect: The LLM's answer is {llm_answer_key}, but the simulated product ({product_canonical_smiles}) does not match any of the provided options."

        except Exception as e:
            return f"An error occurred during the check: {e}"

# Execute the check and print the result.
# result = check_reaction_correctness()
# print(result)