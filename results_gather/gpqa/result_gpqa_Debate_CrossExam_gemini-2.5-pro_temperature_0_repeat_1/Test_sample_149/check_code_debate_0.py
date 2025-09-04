import sys

try:
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("This check requires the RDKit library. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by verifying the reaction pathway.
    """
    # 1. Analyze the starting material based on the LLM's interpretation
    # The LLM correctly identifies the starting material as 4-aminophenylacetaldehyde.
    # Let's verify its molecular formula.
    start_material_smiles = "O=CCc1ccc(N)cc1"  # 4-aminophenylacetaldehyde
    given_formula = "C8H9NO"
    
    mol_start = Chem.MolFromSmiles(start_material_smiles)
    if not mol_start:
        return "Error: Could not parse the SMILES string for the proposed starting material."
        
    calculated_start_formula = CalcMolFormula(mol_start)
    if calculated_start_formula != given_formula:
        return (f"Constraint Violated: The starting material identified by the LLM, "
                f"4-aminophenylacetaldehyde ({calculated_start_formula}), does not match the "
                f"molecular formula given in the question ({given_formula}).")

    # 2. Analyze the intermediate product
    # Reaction 1 & 2: Diazotization + Hydrolysis converts the amine (-NH2) to a hydroxyl (-OH).
    # The intermediate is 4-hydroxyphenylacetaldehyde.
    intermediate_smiles = "O=CCc1ccc(O)cc1"  # 4-hydroxyphenylacetaldehyde
    mol_intermediate = Chem.MolFromSmiles(intermediate_smiles)
    intermediate_formula = CalcMolFormula(mol_intermediate)
    
    # Expected formula: C8H8O2. Let's check.
    if intermediate_formula != "C8H8O2":
        return (f"Error in reasoning: The intermediate product after diazotization/hydrolysis "
                f"should be 4-hydroxyphenylacetaldehyde (C8H8O2), but the derived structure "
                f"has a formula of {intermediate_formula}.")

    # 3. Analyze the final product from Aldol Condensation
    # Reaction 3: Self-aldol condensation (2 molecules react) with heat (dehydration).
    # The overall reaction is: 2 * (C8H8O2) -> Final Product + H2O
    # Expected formula of the final product: C16H16O4 - H2O = C16H14O3.
    expected_final_formula = "C16H14O3"

    # 4. Check the LLM's chosen answer (Option D)
    llm_answer_choice = "D"
    options = {
        "A": "O=CC(c1ccc(O)cc1)C(O)Cc1ccc(O)cc1", # Aldol addition product
        "B": "O=CC(c1ccccc1)C=Cc1ccccc1",
        "C": "O=CCC=Cc1ccc(O)cc1",
        "D": "O=CC(=C(Cc1ccc(O)cc1))c1ccc(O)cc1"  # Aldol condensation product
    }
    
    final_product_smiles = options.get(llm_answer_choice)
    if not final_product_smiles:
        return f"Invalid answer choice '{llm_answer_choice}' provided."

    mol_final = Chem.MolFromSmiles(final_product_smiles)
    calculated_final_formula = CalcMolFormula(mol_final)

    if calculated_final_formula != expected_final_formula:
        return (f"Incorrect Final Product: The chosen answer, Option {llm_answer_choice}, has a molecular "
                f"formula of {calculated_final_formula}. However, the product of a self-aldol "
                f"condensation of 4-hydroxyphenylacetaldehyde should have the formula {expected_final_formula}.")

    # 5. Final check on the logic (Addition vs. Condensation)
    # The LLM correctly states that "Heat" promotes dehydration from the aldol addition product (A)
    # to the aldol condensation product (D). Let's verify the formula of A.
    mol_A = Chem.MolFromSmiles(options["A"])
    formula_A = CalcMolFormula(mol_A)
    if formula_A != "C16H16O4":
        return (f"Error in reasoning: The aldol addition product (Option A) should have the formula "
                f"C16H16O4, but its structure corresponds to {formula_A}.")

    # The LLM's reasoning pathway is chemically and stoichiometrically sound.
    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)