import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def check_answer():
    """
    Checks the correctness of the answer for the reaction of 4,4-dimethylcyclopent-1-enol with bromine.
    """
    try:
        # --- 1. Define Reactants and Products using SMILES for unambiguous representation ---
        # Reactant: 4,4-dimethylcyclopent-1-enol
        reactant_smiles = "CC1(C)CC=C(O)C1"
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        if not reactant_mol:
            return "Error: Could not parse the reactant SMILES string."

        # The provided answer from the LLM
        llm_answer_option = "C"

        # Options provided in the question
        options = {
            "A": "CC1(Br)CCC(=O)C1",  # 4-bromo-4,4-dimethylcyclopentanone
            "B": "CC1(C)CC(Br)C(Br)(O)C1", # (1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol (connectivity)
            "C": "CC1(C)CC(Br)C(=O)C1",  # 2-bromo-4,4-dimethylcyclopentanone
            "D": "CC1(C)CC(Br)C(Br)(O)C1"  # (1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol (connectivity)
        }

        # --- 2. Determine the Expected Product based on Chemical Principles ---
        # Principle: The reaction of an enol with a halogen (Br2) is a classic alpha-halogenation.
        # The reaction proceeds via the enol tautomer. The electron-rich double bond attacks Br2.
        # The bromine adds to the alpha-carbon (the carbon that is not bonded to the -OH group in the enol double bond).
        # The intermediate then loses the proton from the hydroxyl group to reform the very stable carbonyl (C=O) bond.

        # Step-by-step prediction:
        # a. The enol is 4,4-dimethylcyclopent-1-enol. The double bond is between C1 (with OH) and C2.
        # b. The alpha-carbon is C2.
        # c. Bromine will add to C2.
        # d. The C1-OH will become a C1=O (ketone).
        # e. The rest of the structure (4,4-dimethyl group) remains unchanged.
        # f. Expected Product: 2-bromo-4,4-dimethylcyclopentanone.
        expected_product_smiles = "CC1(C)CC(Br)C(=O)C1"
        expected_product_mol = Chem.MolFromSmiles(expected_product_smiles)

        # --- 3. Verify the LLM's Answer ---
        llm_product_smiles = options.get(llm_answer_option)
        if not llm_product_smiles:
            return f"Invalid answer option '{llm_answer_option}'. The options are A, B, C, D."

        llm_product_mol = Chem.MolFromSmiles(llm_product_smiles)

        # Compare the canonical SMILES to check for structural equivalence
        if Chem.MolToSmiles(llm_product_mol, canonical=True) == Chem.MolToSmiles(expected_product_mol, canonical=True):
            # The answer is correct. Let's add a check for why other options are incorrect.
            
            # Check A: 4-bromo-4,4-dimethylcyclopentanone
            mol_A = Chem.MolFromSmiles(options["A"])
            # This would require substitution at a quaternary carbon, which is mechanistically highly unfavorable.
            
            # Check B/D: 1,2-dibromo-4,4-dimethylcyclopentanol
            mol_BD = Chem.MolFromSmiles(options["B"])
            # This is the product of electrophilic addition across the double bond, treating the enol as a simple alkene.
            # While possible, the alpha-halogenation pathway is overwhelmingly favored because it results in the formation
            # of a highly stable carbonyl group, making it the major product.
            
            return "Correct"
        else:
            # The answer is incorrect. Provide the reason.
            reason = (
                f"The provided answer '{llm_answer_option}' corresponds to the product with SMILES '{llm_product_smiles}'.\n"
                f"This is incorrect. The reaction of an enol (4,4-dimethylcyclopent-1-enol) with bromine (Br2) is an alpha-halogenation.\n"
                f"This mechanism involves the addition of a bromine atom to the alpha-carbon (C2) and the subsequent reformation of the stable carbonyl group.\n"
                f"The correct major product is 2-bromo-4,4-dimethylcyclopentanone (SMILES: {expected_product_smiles}), which corresponds to option C."
            )
            # Check if the LLM chose the wrong letter for the right molecule
            if Chem.MolToSmiles(llm_product_mol, canonical=True) != Chem.MolToSmiles(expected_product_mol, canonical=True) and llm_answer_option == "C":
                 return "The LLM chose option C, which is the correct option. However, the reasoning or internal logic might be flawed if it led to a different structure."
            return reason

    except ImportError:
        return "Could not run the check because the 'rdkit' library is not installed. Please install it using 'pip install rdkit-pypi'."
    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check
result = check_answer()
print(result)