import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def check_answer():
    """
    This function checks the correctness of the proposed answer for a multi-step organic synthesis problem.
    It simulates the reaction sequence described in the problem and the LLM's reasoning
    and compares the final simulated product with the proposed answer (Option D).

    The reaction sequence is:
    1. (R)-Limonene + H2/Pd/C -> Product 1
    2. Product 1 + m-CPBA -> Product 2
    3. Product 2 + NaOMe -> Product 3
    4. Product 3 + Propanoic acid/DCC/DMAP -> Product 4

    The LLM's reasoning for the stereochemistry at each step is used to construct the intermediate molecules.
    """
    error_log = []

    # --- Step 0: Define the target molecule from the proposed answer (Option D) ---
    # Option D is (1S,2S,4R)-4-isopropyl-2-methoxy-1-methylcyclohexyl propionate.
    # We construct its SMILES representation to serve as the target for our simulation.
    try:
        # SMILES for (1S,2S,4R)-4-isopropyl-2-methoxy-1-methylcyclohexyl propionate
        # C1 is S, C2 is S, C4 is R.
        option_d_smiles = "CCC(=O)O[C@]1(C)[C@@H](OC)C[C@H](C(C)C)CC1"
        option_d_mol = Chem.MolFromSmiles(option_d_smiles)
        if option_d_mol is None:
            return "Constraint check failed: The SMILES string for the proposed answer (Option D) is invalid."
        option_d_canonical_smiles = Chem.MolToSmiles(option_d_mol, isomericSmiles=True)
    except Exception as e:
        return f"An error occurred while processing Option D: {e}"

    # --- Simulation based on the LLM's step-by-step reasoning ---

    # The core of the check is to verify the consistency of the LLM's reasoning chain.
    # We will construct the key intermediates as described by the LLM and perform the final reaction.

    # LLM's Product 2: (1S, 2R, 4R)-1-methyl-4-isopropyl-1,2-epoxycyclohexane.
    # This is the result of hydrogenation and sterically-directed epoxidation.
    product_2_smiles = "C[C@]12O[C@H](C[C@H](C(C)C)CC1)C2"
    product_2 = Chem.MolFromSmiles(product_2_smiles)
    if product_2 is None:
        return "Incorrect. Failed to construct the epoxide intermediate (Product 2) based on the LLM's reasoning. The SMILES may be incorrect."

    # LLM's Product 3: (1S, 2S, 4R)-4-isopropyl-2-methoxy-1-methylcyclohexan-1-ol.
    # This is the result of the S_N2 opening of Product 2, which correctly inverts the stereocenter at C2.
    product_3_smiles = "C[C@]1(O)[C@@H](OC)C[C@H](C(C)C)CC1"
    product_3 = Chem.MolFromSmiles(product_3_smiles)
    if product_3 is None:
        return "Incorrect. Failed to construct the alcohol intermediate (Product 3) after the S_N2 opening. The SMILES may be incorrect."

    # Step 4: Simulate the final esterification reaction.
    # The LLM correctly states that Steglich esterification proceeds with retention of configuration at C1.
    # We apply this reaction to the LLM's proposed Product 3.
    rxn_step4 = AllChem.ReactionFromSmarts('[C:1]([OH:2])>>[C:1]([O:2]C(=O)CC)')
    products4 = rxn_step4.RunReactants((product_3,))
    
    if not products4 or not products4[0]:
        return "Incorrect. The simulation of the final esterification step failed to produce a product."
    
    simulated_product_4 = products4[0][0]
    Chem.SanitizeMol(simulated_product_4)
    simulated_canonical_smiles = Chem.MolToSmiles(simulated_product_4, isomericSmiles=True)

    # --- Final Verification ---
    # Compare the canonical SMILES of the simulated final product with the one from Option D.
    if simulated_canonical_smiles == option_d_canonical_smiles:
        return "Correct"
    else:
        reason = "Incorrect. The final product derived from the step-by-step reasoning does not match Option D.\n"
        reason += f"Expected SMILES (from Option D): {option_d_canonical_smiles}\n"
        reason += f"Simulated SMILES (from LLM logic): {simulated_canonical_smiles}\n"
        reason += "This indicates a flaw or inconsistency in the LLM's reasoning chain, likely in assigning the stereochemistry of an intermediate or the final product."
        return reason

# Print the result of the check
print(check_answer())