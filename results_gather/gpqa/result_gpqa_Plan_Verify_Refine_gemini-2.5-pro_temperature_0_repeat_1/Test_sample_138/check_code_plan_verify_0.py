import sys
import io

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer for the given chemistry problem.
    It verifies the transformation of reactants to products based on the principles
    of alpha-oxidation of ketones.
    """
    # Suppress RDKit's non-critical warnings for cleaner output
    original_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        # Attempt to import necessary libraries
        from rdkit import Chem
        import pubchempy as pcp
    except ImportError:
        # If libraries are not found, return an informative error message
        return ("Required libraries not found. Please install rdkit-pypi and pubchempy "
                "by running: pip install rdkit-pypi pubchempy")
    finally:
        # Restore standard error stream
        sys.stderr = original_stderr

    # --- Problem and Answer Definition ---
    llm_answer = "C"

    options = {
        "A": {"A": "4-isopropylcyclohexan-1-one", "B": "5-methylhexane-2,3-diol"},
        "B": {"A": "4-isopropyl-2-methoxycyclohexan-1-ol", "B": "5-methylhexane-2,3-diol"},
        "C": {"A": "4-isopropylcyclohexan-1-one", "B": "5-methylhexan-2-one"},
        "D": {"A": "4-isopropyl-2-methoxycyclohexan-1-ol", "B": "5-methylhexan-2-one"}
    }

    target_products = {
        "A": "4-isopropylcyclohexane-1,2-dione",
        "B": "5-methylhexane-2,3-dione"
    }

    # --- Helper Functions ---
    memo = {}  # Cache for PubChem lookups to improve performance
    def get_mol_from_name(name):
        """Retrieves a canonical RDKit molecule object from its chemical name using PubChem."""
        if name in memo:
            return memo[name]
        try:
            compounds = pcp.get_compounds(name, 'name')
            if not compounds:
                result = (None, f"Could not find compound '{name}' in PubChem.")
            else:
                smiles = compounds[0].canonical_smiles
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    result = (None, f"Could not parse SMILES for '{name}'.")
                else:
                    result = (mol, None)
        except Exception as e:
            result = (None, f"Network or API error fetching '{name}': {e}")
        
        memo[name] = result
        return result

    def check_transformation(reactant_name, product_name):
        """
        Checks if the reactant can form the product via alpha-oxidation of a ketone.
        This is done by checking the reverse reaction: reducing one of the alpha-carbonyls
        of the product to a methylene group and seeing if it matches the reactant.
        """
        # 1. Rule-based check for invalid reactant types for this specific reaction
        if 'diol' in reactant_name.lower() or 'ol' in reactant_name.lower():
            return False, f"Reactant '{reactant_name}' is an alcohol/diol. The reaction (alpha-nitrosation) requires a ketone with an alpha-methylene or alpha-methine group."

        # 2. Get molecule structures from names
        reactant_mol, error = get_mol_from_name(reactant_name)
        if error: return False, error
        product_mol, error = get_mol_from_name(product_name)
        if error: return False, error

        reactant_smiles = Chem.MolToSmiles(reactant_mol, canonical=True)

        # 3. Find alpha-diketone groups in the product to identify reaction centers
        diketone_core_pattern = Chem.MolFromSmarts('[C:1](=[O])-[C:2](=[O])')
        diketone_matches = product_mol.GetSubstructMatches(diketone_core_pattern)

        if not diketone_matches:
            return False, f"Product '{product_name}' is not an alpha-diketone as expected."

        # 4. Generate possible precursors by simulating the reverse reaction (reduction)
        possible_precursors_smiles = set()
        for c1_idx, c2_idx in diketone_matches:
            # Try reducing the first carbonyl (C1) to a methylene (CH2)
            mol1 = Chem.RWMol(product_mol)
            atom1 = mol1.GetAtomWithIdx(c1_idx)
            atom1.SetAtomicNum(6)
            atom1.SetNumExplicitHs(2)
            atom1.SetNoImplicit(True)
            for neighbor in list(atom1.GetNeighbors()):
                if neighbor.GetAtomicNum() == 8 and mol1.GetBondBetweenAtoms(c1_idx, neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                    mol1.RemoveBond(c1_idx, neighbor.GetIdx())
                    mol1.RemoveAtom(neighbor.GetIdx())
                    break
            try:
                Chem.SanitizeMol(mol1)
                possible_precursors_smiles.add(Chem.MolToSmiles(Chem.RemoveHs(mol1), canonical=True))
            except Chem.rdchem.MolSanitizeException: pass

            # Try reducing the second carbonyl (C2) to a methylene (CH2)
            mol2 = Chem.RWMol(product_mol)
            atom2 = mol2.GetAtomWithIdx(c2_idx)
            atom2.SetAtomicNum(6)
            atom2.SetNumExplicitHs(2)
            atom2.SetNoImplicit(True)
            for neighbor in list(atom2.GetNeighbors()):
                if neighbor.GetAtomicNum() == 8 and mol2.GetBondBetweenAtoms(c2_idx, neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                    mol2.RemoveBond(c2_idx, neighbor.GetIdx())
                    mol2.RemoveAtom(neighbor.GetIdx())
                    break
            try:
                Chem.SanitizeMol(mol2)
                possible_precursors_smiles.add(Chem.MolToSmiles(Chem.RemoveHs(mol2), canonical=True))
            except Chem.rdchem.MolSanitizeException: pass

        # 5. Check if the given reactant is among the valid, generated precursors
        if reactant_smiles in possible_precursors_smiles:
            return True, ""
        else:
            return False, f"The structure of '{reactant_name}' does not lead to '{product_name}' via alpha-oxidation. Working backwards, reducing the product's carbonyls does not yield the reactant."

    # --- Main Verification Logic ---
    correct_option = None
    failure_reasons = {}

    for option_key, reactants in options.items():
        reactant_A_name = reactants["A"]
        reactant_B_name = reactants["B"]
        
        is_A_correct, reason_A = check_transformation(reactant_A_name, target_products["A"])
        is_B_correct, reason_B = check_transformation(reactant_B_name, target_products["B"])

        if is_A_correct and is_B_correct:
            correct_option = option_key
            break  # Found the correct option, no need to check further
        else:
            # Store reasons for failure for detailed feedback
            reasons = []
            if not is_A_correct: reasons.append(f"For reactant A ({reactant_A_name}): {reason_A}")
            if not is_B_correct: reasons.append(f"For reactant B ({reactant_B_name}): {reason_B}")
            failure_reasons[option_key] = " ".join(reasons)

    # --- Final Verdict ---
    if correct_option == llm_answer:
        return "Correct"
    elif correct_option is None:
        reason = failure_reasons.get(llm_answer, "The code could not validate any option.")
        return f"Incorrect. The LLM's answer is {llm_answer}, but the code could not find any valid option. The reason the LLM's choice is wrong is: {reason}"
    else:
        llm_reason = failure_reasons.get(llm_answer, "The proposed reactants do not yield the target products.")
        return (f"Incorrect. The LLM's answer is {llm_answer}, but this is wrong. "
                f"Reason: {llm_reason} "
                f"The correct option is {correct_option} because both of its proposed reactants correctly yield the target products via alpha-oxidation.")

# Execute the check and print the result to the console.
result = check_chemistry_answer()
print(result)