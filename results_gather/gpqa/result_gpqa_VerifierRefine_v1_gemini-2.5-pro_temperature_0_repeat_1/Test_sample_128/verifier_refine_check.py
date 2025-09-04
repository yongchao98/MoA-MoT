def check_correctness():
    """
    Checks the correctness of the LLM's answer for the multi-step reaction problem.
    This function uses the RDKit library to model and verify each chemical step.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
    except ImportError:
        return "Could not check the answer because the 'rdkit' library is not installed. Please install it using 'pip install rdkit'."

    errors = []
    llm_answer_text = "B) 3,4-dimethylcyclohexan-1-one"
    llm_final_choice = "B"

    # --- Verification Steps ---

    # Step 1: Verify Compound A from Hint (a)
    # Hint (a) implies a Wittig reaction. The LLM correctly performs a retrosynthesis.
    # Product: 1,2-dimethyl-4-(propan-2-ylidene)cyclopentane
    # LLM's proposed Compound A: 3,4-dimethylcyclopentan-1-one
    # The LLM's interpretation is chemically sound because naming the ketone gives the carbonyl C1 priority,
    # leading to 3,4-dimethylcyclopentan-1-one. This logic is correct.
    try:
        compound_A_smiles = "CC1CC(=O)C(C)C1"  # 3,4-dimethylcyclopentan-1-one
        compound_A_mol = Chem.MolFromSmiles(compound_A_smiles)
        if compound_A_mol is None: raise ValueError("Invalid SMILES for Compound A")
    except Exception as e:
        errors.append(f"Failed to validate Step 1 (Compound A derivation): {e}")

    # Step 2: Verify Compound A with IR Hint (b)
    # Hint (b) says A has a peak at ~1750 cm-1, characteristic of a cyclopentanone (5-membered ring ketone).
    try:
        compound_A_mol = Chem.MolFromSmiles("CC1CC(=O)C(C)C1")
        carbonyl_atom = next((atom for atom in compound_A_mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1), None)
        if carbonyl_atom is None:
            errors.append("Proposed Compound A does not contain a carbonyl group.")
        else:
            carbonyl_carbon = carbonyl_atom.GetNeighbors()[0]
            if not carbonyl_carbon.IsInRingSize(5):
                errors.append("Constraint check failed for Compound A: The IR hint (~1750 cm-1) implies a 5-membered ring ketone, but the proposed structure is not a cyclopentanone.")
    except Exception as e:
        errors.append(f"Failed to validate Step 2 (IR hint for A): {e}")

    # Step 3: Verify the reaction sequence A -> B -> C -> E
    try:
        # A -> B (HCN addition)
        llm_B_smiles = "CC1CC(O)(C#N)C(C)C1"  # 1-cyano-3,4-dimethylcyclopentan-1-ol
        rxn_A_to_B = AllChem.ReactionFromSmarts('[C:1]=[O:2]>>[C:1](O)[C#N]')
        products = rxn_A_to_B.RunReactants((Chem.MolFromSmiles("CC1CC(=O)C(C)C1"),))
        if not products or Chem.MolToSmiles(products[0][0], True) != Chem.MolToSmiles(Chem.MolFromSmiles(llm_B_smiles), True):
            errors.append("The structure of Compound B derived from Compound A is incorrect.")

        # B -> C (Nitrile reduction)
        llm_C_smiles = "CC1CC(O)(CN)C(C)C1"  # 1-(aminomethyl)-3,4-dimethylcyclopentan-1-ol
        rxn_B_to_C = AllChem.ReactionFromSmarts('[C:1]#[N:2]>>[C:1]C[N:2]')
        products = rxn_B_to_C.RunReactants((Chem.MolFromSmiles(llm_B_smiles),))
        if not products or Chem.MolToSmiles(products[0][0], True) != Chem.MolToSmiles(Chem.MolFromSmiles(llm_C_smiles), True):
            errors.append("The structure of Compound C derived from Compound B is incorrect.")

        # C -> E (Tiffeneau-Demjanov rearrangement)
        # This involves ring expansion from a 5-membered ring to a 6-membered ring.
        compound_C_mol = Chem.MolFromSmiles(llm_C_smiles)
        llm_E_smiles = "CC1CCC(=O)C(C)C1"  # 3,4-dimethylcyclohexan-1-one
        compound_E_mol = Chem.MolFromSmiles(llm_E_smiles)
        
        is_c_cyclopentanol_deriv = any(atom.IsInRingSize(5) for atom in compound_C_mol.GetAtoms())
        is_e_cyclohexanone_deriv = any(atom.IsInRingSize(6) for atom in compound_E_mol.GetAtoms())

        if not (is_c_cyclopentanol_deriv and is_e_cyclohexanone_deriv):
            errors.append("The C -> E transformation did not result in the expected ring expansion from a 5-membered to a 6-membered ring.")
        
        # Check molecular formulas (C8H17NO -> C8H14O)
        formula_C = rdMolDescriptors.CalcMolFormula(compound_C_mol)
        formula_E = rdMolDescriptors.CalcMolFormula(compound_E_mol)
        if formula_C != "C8H17NO" or formula_E != "C8H14O":
            errors.append(f"Molecular formula mismatch during C -> E step. C should be C8H17NO (is {formula_C}), E should be C8H14O (is {formula_E}).")

    except Exception as e:
        errors.append(f"Failed to validate Step 3 (Reaction sequence): {e}")

    # Step 4: Verify Compound E with IR Hint (b)
    # Hint (b) says E has a peak at ~1715 cm-1, characteristic of a cyclohexanone (6-membered ring ketone).
    try:
        compound_E_mol = Chem.MolFromSmiles("CC1CCC(=O)C(C)C1")
        carbonyl_atom = next((atom for atom in compound_E_mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1), None)
        if carbonyl_atom is None:
            errors.append("Proposed Compound E does not contain a carbonyl group.")
        else:
            carbonyl_carbon = carbonyl_atom.GetNeighbors()[0]
            if not carbonyl_carbon.IsInRingSize(6):
                errors.append("Constraint check failed for Compound E: The IR hint (~1715 cm-1) implies a 6-membered ring ketone, but the proposed structure is not a cyclohexanone.")
    except Exception as e:
        errors.append(f"Failed to validate Step 4 (IR hint for E): {e}")

    # Final check: Does the derived structure match the selected option?
    if "3,4-dimethylcyclohexan-1-one" not in llm_answer_text:
        errors.append(f"The final answer text '{llm_answer_text}' does not match the derived structure '3,4-dimethylcyclohexan-1-one'.")
    
    if llm_final_choice != "B":
        errors.append(f"The final choice <<<B>>> is incorrect. The correct option corresponding to 3,4-dimethylcyclohexan-1-one is B.")


    if not errors:
        return "Correct"
    else:
        return "The answer is incorrect. Here are the reasons:\n- " + "\n- ".join(errors)

# Execute the check and print the result
result = check_correctness()
print(result)