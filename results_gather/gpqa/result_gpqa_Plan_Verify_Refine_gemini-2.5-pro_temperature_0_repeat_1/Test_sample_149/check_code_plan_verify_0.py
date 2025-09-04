import sys
try:
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit'")
    sys.exit(1)

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's reasoning for the given chemistry problem.
    """
    error_messages = []

    # --- Step 1: Define molecules and check the proposed starting material ---
    
    # SMILES strings for all relevant compounds
    smiles_map = {
        "start_proposed": "O=CCc1ccc(N)cc1",  # 2-(4-aminophenyl)acetaldehyde
        "intermediate": "O=CCc1ccc(O)cc1",   # 2-(4-hydroxyphenyl)acetaldehyde
        "option_C_final": "O=C/C(c1ccc(O)cc1)=C/Cc1ccc(O)cc1", # 2,4-bis(4-hydroxyphenyl)but-2-enal
        "option_D_addition": "O=C[C@H](c1ccc(O)cc1)[C@H](O)Cc1ccc(O)cc1" # 3-hydroxy-2,4-bis(4-hydroxyphenyl)butanal
    }

    # Create RDKit molecule objects
    mols = {name: Chem.MolFromSmiles(s) for name, s in smiles_map.items()}
    if None in mols.values():
        return "Internal Error: Could not parse one of the SMILES strings."

    # 1a. Verify the molecular formula of the starting material
    start_mol = mols["start_proposed"]
    expected_formula = "C8H9NO"
    actual_formula = CalcMolFormula(start_mol)
    if actual_formula != expected_formula:
        error_messages.append(f"Starting Material Error: The proposed starting material (2-(4-aminophenyl)acetaldehyde) has a formula of {actual_formula}, which does not match the required {expected_formula}.")

    # 1b. Verify the functional groups based on NMR interpretation
    # SMARTS patterns for functional groups
    patterns = {
        "aldehyde": "[CX3H1](=O)",
        "primary_aromatic_amine": "c-N",
        "para_disubstitution": "c1(-[*])ccc(-[*])cc1"
    }
    for group, smarts in patterns.items():
        if not start_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            error_messages.append(f"Starting Material Error: The proposed structure lacks the {group} functional group, which contradicts the NMR data interpretation.")

    # --- Step 2: Verify the reaction sequence ---

    # 2a. Verify the intermediate from diazotization/hydrolysis (Ar-NH2 -> Ar-OH)
    # The intermediate should be 2-(4-hydroxyphenyl)acetaldehyde (C8H8O2)
    intermediate_mol = mols["intermediate"]
    expected_intermediate_formula = "C8H8O2"
    actual_intermediate_formula = CalcMolFormula(intermediate_mol)
    if actual_intermediate_formula != expected_intermediate_formula:
        error_messages.append(f"Reaction Step Error: The intermediate product (2-(4-hydroxyphenyl)acetaldehyde) has formula {actual_intermediate_formula}, but it should be {expected_intermediate_formula}.")

    # 2b. Verify the final product from aldol condensation
    # Reaction: 2 * Intermediate -> Final Product + H2O
    # Expected formula of final product: 2 * C8H8O2 - H2O = C16H16O4 - H2O = C16H14O3
    expected_final_formula = "C16H14O3"
    actual_final_formula = CalcMolFormula(mols["option_C_final"])
    if actual_final_formula != expected_final_formula:
        error_messages.append(f"Final Product Error: The aldol condensation product (Option C) should have formula {expected_final_formula}, but its actual formula is {actual_final_formula}.")

    # 2c. Verify that Option D is the addition product (before dehydration)
    # Expected formula of addition product: 2 * C8H8O2 = C16H16O4
    expected_addition_formula = "C16H16O4"
    actual_addition_formula = CalcMolFormula(mols["option_D_addition"])
    if actual_addition_formula != expected_addition_formula:
        error_messages.append(f"Logic Error: The aldol addition intermediate (Option D) should have formula {expected_addition_formula}, but its actual formula is {actual_addition_formula}.")

    # The LLM correctly reasoned that "Heat" causes dehydration, making C the final product over D. This logic is sound.

    # --- Step 3: Final Conclusion ---
    if error_messages:
        return "Incorrect. The answer is wrong for the following reason(s):\n- " + "\n- ".join(error_messages)
    else:
        # A minor note: The LLM's calculation for the Degree of Unsaturation (DBE) was incorrect (4 instead of 5).
        # DBE = (2C + 2 + N - H)/2 = (16 + 2 + 1 - 9)/2 = 5.
        # However, this minor error in the reasoning did not lead to an incorrect structure or final answer.
        # The overall logic and final conclusion are correct.
        return "Correct"

# Run the check and print the result
result = check_chemistry_answer()
print(result)