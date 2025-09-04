try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit'")
    # Create a dummy function to avoid crashing if rdkit is not installed
    def run_check():
        return "Skipped: RDKit library is required for this check but was not found."
else:
    def get_mol_formula(mol):
        """Calculates the molecular formula of an RDKit molecule object."""
        if mol is None:
            return None
        return Chem.rdMolDescriptors.CalcMolFormula(mol)

    def check_substructure(mol, smarts):
        """Checks if a molecule contains a given substructure (defined by SMARTS)."""
        if mol is None:
            return False
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            # This indicates an error in the SMARTS pattern itself
            return False
        return mol.HasSubstructMatch(pattern)

    def run_check():
        """
        Checks the correctness of the answer by verifying the chemical logic.
        """
        # --- Define Structures and Constraints ---
        
        # Constraint from the question
        target_formula = "C11H12O"
        
        # Product identified from NMR data in the reasoning: (E)-4-(p-tolyl)but-3-en-2-one
        product_smiles = "CC(=O)/C=C/c1ccc(C)cc1"
        product_mol = Chem.MolFromSmiles(product_smiles)

        # Options for Compound X as given in the question
        options = {
            "A": {"name": "2-(4-methylstyryl)oxirane", "smiles": "Cc1ccc(/C=C/C2OC2)cc1"},
            "B": {"name": "2-styrylepoxide", "smiles": "c1ccc(/C=C/C2OC2)cc1"},
            "C": {"name": "2-methyl-3-styryloxirane", "smiles": "CC1(OC1/C=C/c2ccccc2)"},
            "D": {"name": "2-(1-phenylprop-1-en-2-yl)oxirane", "smiles": "CC(=C/c1ccccc1)C2OC2"}
        }

        # SMARTS patterns for key functional groups identified in the reasoning
        smarts_patterns = {
            "p_tolyl": "c1(C)ccc(cc1)-*",
            "ketone": "[#6][CX3](=O)[#6]",
            "epoxide": "C1OC1"
        }
        
        # --- Step 1: Verify the identified product structure ---
        # The reasoning's first step is to identify the product. Let's check if the proposed
        # product structure is consistent with the key features from the NMR data.
        
        product_formula = get_mol_formula(product_mol)
        if product_formula != target_formula:
            return f"Incorrect. The reasoning is flawed because the identified product ({product_formula}) does not have the required molecular formula {target_formula}."

        if not check_substructure(product_mol, smarts_patterns["ketone"]):
            return "Incorrect. The reasoning is flawed because the identified product does not contain a ketone, which contradicts the 13C NMR signal at 197.7 ppm."

        if not check_substructure(product_mol, smarts_patterns["p_tolyl"]):
            return "Incorrect. The reasoning is flawed because the identified product does not contain a p-tolyl group, which contradicts the 1H NMR data."

        # --- Step 2: Evaluate all options against the logical constraints ---
        # The core of the reasoning is that the starting material must have the same formula
        # and the same p-tolyl skeleton as the product.
        
        valid_candidates = []
        for option_key, data in options.items():
            mol = Chem.MolFromSmiles(data["smiles"])
            if mol is None:
                continue # Skip invalid SMILES

            # Check if the option meets all criteria
            formula_ok = (get_mol_formula(mol) == target_formula)
            skeleton_ok = check_substructure(mol, smarts_patterns["p_tolyl"])
            is_epoxide = check_substructure(mol, smarts_patterns["epoxide"])

            if formula_ok and skeleton_ok and is_epoxide:
                valid_candidates.append(option_key)

        # --- Step 3: Compare the result with the provided answer ---
        final_answer = "A"
        
        if final_answer in valid_candidates and len(valid_candidates) == 1:
            return "Correct"
        elif final_answer not in valid_candidates:
            # Find out why the proposed answer 'A' failed the check
            mol_A = Chem.MolFromSmiles(options["A"]["smiles"])
            reasons = []
            if get_mol_formula(mol_A) != target_formula:
                reasons.append(f"its formula is {get_mol_formula(mol_A)}, not {target_formula}")
            if not check_substructure(mol_A, smarts_patterns["p_tolyl"]):
                reasons.append("it does not contain the required p-tolyl group")
            if not check_substructure(mol_A, smarts_patterns["epoxide"]):
                reasons.append("it is not an epoxide")
            
            return f"Incorrect. The provided answer '{final_answer}' is wrong because " + ", ".join(reasons) + "."
        else: # More than one candidate is valid
            return f"Incorrect. The reasoning is not selective enough, as options {valid_candidates} all fit the criteria. The provided answer '{final_answer}' may be correct, but the logic is incomplete."

# Run the check and print the result
print(run_check())