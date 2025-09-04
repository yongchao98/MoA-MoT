import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def check_answer():
    """
    Checks the correctness of the provided answer by verifying molecular formulas
    and structural features of the compounds involved.
    """
    # --- Define Molecules and Constraints from the Question ---

    # Constraint 1: Molecular formula of the starting material (Compound X)
    reactant_formula_str = "C11H12O"

    # The product structure deduced from NMR data in the answer
    # Product: 4-(4-methylphenyl)but-3-en-2-one
    product_smiles = "CC(=O)C=Cc1ccc(C)cc1"
    product_mol = Chem.MolFromSmiles(product_smiles)
    if not product_mol:
        return "Error: Could not parse the product SMILES string."

    # The options for Compound X
    options = {
        "A": {"name": "2-styrylepoxide", "smiles": "c1ccc(cc1)C=CC2OC2"},
        "B": {"name": "2-(1-phenylprop-1-en-2-yl)oxirane", "smiles": "CC(=Cc1ccccc1)C2OC2"},
        "C": {"name": "2-methyl-3-styryloxirane", "smiles": "CC1OC1C=Cc1ccccc1"},
        "D": {"name": "2-(4-methylstyryl)oxirane", "smiles": "Cc1ccc(C=CC2OC2)cc1"}
    }

    # --- Verification Steps ---

    # 1. Verify the properties of the deduced product
    # The reaction is an isomerization, so the product should have the same formula as the reactant.
    product_formula = CalcMolFormula(product_mol)
    if product_formula != reactant_formula_str:
        return (f"Incorrect: The answer's deduced product ({product_formula}) does not have the "
                f"same molecular formula as the starting material ({reactant_formula_str}). "
                f"This contradicts the assumption of an isomerization reaction.")

    # 2. Define substructures for checking
    # p-tolyl group: a benzene ring with a methyl group, para to the rest of the molecule.
    # We can check for a methyl-substituted benzene ring.
    p_tolyl_pattern = Chem.MolFromSmarts("c1ccc(C)cc1")
    # Phenyl group: an unsubstituted benzene ring.
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")

    # 3. Check if the product contains the p-tolyl group, as reasoned in the answer.
    if not product_mol.HasSubstructMatch(p_tolyl_pattern):
        return ("Incorrect: The answer's deduced product (4-(4-methylphenyl)but-3-en-2-one) "
                "does not contain a p-tolyl group, which contradicts the NMR analysis.")

    # 4. Evaluate each option for Compound X
    correct_option = None
    errors = []

    for key, data in options.items():
        mol = Chem.MolFromSmiles(data["smiles"])
        if not mol:
            errors.append(f"Option {key}: Could not parse SMILES '{data['smiles']}'.")
            continue

        # Check 1: Molecular Formula
        formula = CalcMolFormula(mol)
        if formula != reactant_formula_str:
            if key == "D": # The proposed correct answer
                 return f"Incorrect: The answer claims D is correct, but its formula ({formula}) does not match the required {reactant_formula_str}."
            continue # This option is invalid, move to the next one.

        # Check 2: Presence of the required p-tolyl group
        # Since the product has a p-tolyl group and the reaction is an isomerization,
        # the reactant must also have this group.
        has_p_tolyl = mol.HasSubstructMatch(p_tolyl_pattern)
        
        # To be precise, we check that it's not just a phenyl group.
        # A p-tolyl group will match both patterns, but a phenyl group will only match the phenyl pattern.
        has_phenyl_only = mol.HasSubstructMatch(phenyl_pattern) and not has_p_tolyl

        if has_p_tolyl:
            if correct_option is None:
                correct_option = key
            else:
                # This would mean more than one option is plausible
                return f"Incorrect: The reasoning is flawed. Both option {correct_option} and {key} have the correct formula and a p-tolyl group."
        
        elif has_phenyl_only:
            if key == "D": # The proposed correct answer
                return f"Incorrect: The answer claims D is correct, but it only contains a phenyl group, not the required p-tolyl group found in the product."
        else:
            # This case shouldn't happen for the given options but is good for robustness
            if key == "D":
                return f"Incorrect: The answer claims D is correct, but it contains neither a phenyl nor a p-tolyl group."


    # --- Final Conclusion ---
    if correct_option == "D":
        # Verify the answer's specific claims about other options
        # Claim: A has wrong formula.
        mol_A = Chem.MolFromSmiles(options["A"]["smiles"])
        if CalcMolFormula(mol_A) == reactant_formula_str:
            return "Incorrect: The answer's reasoning that option A has the wrong formula is false."
        
        # Claim: B and C have phenyl groups, not p-tolyl.
        mol_B = Chem.MolFromSmiles(options["B"]["smiles"])
        mol_C = Chem.MolFromSmiles(options["C"]["smiles"])
        if not (mol_B.HasSubstructMatch(phenyl_pattern) and not mol_B.HasSubstructMatch(p_tolyl_pattern)):
             return "Incorrect: The answer's reasoning that option B has a phenyl group (and not p-tolyl) is false."
        if not (mol_C.HasSubstructMatch(phenyl_pattern) and not mol_C.HasSubstructMatch(p_tolyl_pattern)):
             return "Incorrect: The answer's reasoning that option C has a phenyl group (and not p-tolyl) is false."

        return "Correct"
    elif correct_option is None:
        return "Incorrect: No option satisfies all the constraints (correct formula and p-tolyl group)."
    else:
        return f"Incorrect: The code identified option {correct_option} as correct, but the provided answer was D."

# Run the check and print the result
result = check_answer()
print(result)