import sys
import importlib

def check_chemistry_answer():
    """
    This function checks the correctness of the provided multi-step chemistry problem.
    It verifies each step of the reaction sequence and then uses the RDKit library
    to analyze the final product's structure and count its chemically distinct hydrogens.
    """
    # --- Step 0: Check for RDKit installation ---
    try:
        importlib.import_module("rdkit")
    except ImportError:
        return ("The 'rdkit' library is required to run this check but is not installed.\n"
                "Please install it, for example, using pip: pip install rdkit")

    from rdkit import Chem

    # --- Step 1: Define the expected reaction pathway and final answer from the LLM ---
    # The LLM's answer identifies the products and the final count.
    llm_reasoning = {
        "product_1": "2-bromocyclohexanone",
        "product_2": "Cyclopentanecarboxylic acid",
        "product_3": "Cyclopentanecarbonyl chloride",
        "product_4": "Cyclopentanecarbaldehyde",
        "final_count": 6
    }

    # --- Step 2: Verify the logic of the reaction pathway ---
    # The code will check if the LLM's identification of each product is correct based on standard organic chemistry rules.
    
    # Reaction 1: Alpha-bromination of a ketone. Correctly identified.
    if llm_reasoning["product_1"] != "2-bromocyclohexanone":
        return "Incorrect identification of Product 1. The alpha-bromination of cyclohexanone yields 2-bromocyclohexanone."

    # Reaction 2: Favorskii rearrangement. This is a key step involving ring contraction. Correctly identified.
    if llm_reasoning["product_2"] != "Cyclopentanecarboxylic acid":
        return "Incorrect identification of Product 2. The reaction is a Favorskii rearrangement, which contracts the 6-membered ring to a 5-membered ring, yielding Cyclopentanecarboxylic acid."

    # Reaction 3: Conversion of carboxylic acid to acyl chloride. Correctly identified.
    if llm_reasoning["product_3"] != "Cyclopentanecarbonyl chloride":
        return "Incorrect identification of Product 3. The reaction of a carboxylic acid with thionyl chloride yields an acyl chloride, in this case, Cyclopentanecarbonyl chloride."

    # Reaction 4: Selective reduction of acyl chloride to aldehyde. Correctly identified.
    if llm_reasoning["product_4"] != "Cyclopentanecarbaldehyde":
        return "Incorrect identification of Product 4. Lithium tri-tert-butoxyaluminum hydride is a selective reducing agent that converts acyl chlorides to aldehydes, yielding Cyclopentanecarbaldehyde."

    # --- Step 3: Analyze the final product using RDKit ---
    # If the reaction pathway is correct, we now verify the count of distinct hydrogens on the final product.
    
    # SMILES string for the final product, Cyclopentanecarbaldehyde
    product_4_smiles = "O=CC1CCCC1"
    
    try:
        mol = Chem.MolFromSmiles(product_4_smiles)
        if mol is None:
            return f"Error: RDKit could not parse the SMILES for the final product: {product_4_smiles}"

        # Add explicit hydrogens to the molecule graph
        mol_with_hs = Chem.AddHs(mol)

        # Use canonical atom ranking to find symmetrically equivalent atoms.
        # Atoms with the same rank are considered chemically equivalent.
        ranks = list(Chem.CanonicalRankAtoms(mol_with_hs, breakTies=False))

        # Collect the ranks of all hydrogen atoms
        hydrogen_ranks = set()
        for atom in mol_with_hs.GetAtoms():
            if atom.GetAtomicNum() == 1:  # Atomic number for Hydrogen is 1
                hydrogen_ranks.add(ranks[atom.GetIdx()])

        # The number of distinct hydrogens is the number of unique ranks found.
        calculated_distinct_hydrogens = len(hydrogen_ranks)

        # Compare the calculated count with the LLM's answer
        if calculated_distinct_hydrogens == llm_reasoning["final_count"]:
            return "Correct"
        else:
            return (f"The final count of distinct hydrogens is incorrect. "
                    f"The code calculated {calculated_distinct_hydrogens} distinct hydrogen environments for Cyclopentanecarbaldehyde, "
                    f"while the provided answer states {llm_reasoning['final_count']}.")

    except Exception as e:
        return f"An unexpected error occurred during the RDKit analysis: {e}"

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)