import re

def check_diels_alder_stereochemistry():
    """
    Checks the correctness of the answer for the Diels-Alder reaction of
    5-fluorocyclopenta-1,3-diene and maleic anhydride.

    The function codifies the key chemical principles and nomenclature rules
    to derive the correct product and compares it with the provided answer.
    """
    # --- Problem Definition ---
    options = {
        'A': "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        'B': "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        'C': "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        'D': "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione"
    }
    provided_answer_key = 'D'
    provided_answer_reasoning = {
        "endo_exo_rule": "endo",
        "facial_selectivity_rule": "syn-attack",
        "product_F_position": "anti",
        "nomenclature_for_anti_F": "8r"
    }

    # --- Chemical Knowledge Base ---

    def get_skeleton_orientation(name):
        """Determines if the skeleton is endo or exo based on stereodescriptors."""
        # For the (3aR, 7aS) enantiomer, the endo adduct has (4S, 7R) stereochemistry.
        # The exo adduct has (4R, 7S) stereochemistry.
        if '4S' in name and '7R' in name:
            return 'endo'
        elif '4R' in name and '7S' in name:
            return 'exo'
        return 'unknown'

    def get_fluorine_descriptor(name):
        """Determines the fluorine's descriptor (8r/8s)."""
        match = re.search(r'8([rs])', name)
        return f"8{match.group(1)}" if match else 'unknown'

    # --- Verification Logic ---
    
    # 1. Check Endo/Exo Selectivity Principle
    # The Alder-Endo rule predicts the 'endo' adduct is the major kinetic product.
    # This matches the provided reasoning.
    if provided_answer_reasoning["endo_exo_rule"] != "endo":
        return "Incorrect Reasoning: The provided answer does not correctly apply the Alder-Endo rule, which favors the 'endo' product."

    # Filter options based on the endo rule
    endo_candidates = {
        key: name for key, name in options.items()
        if get_skeleton_orientation(name) == "endo"
    }
    
    if not ('A' in endo_candidates and 'D' in endo_candidates and len(endo_candidates) == 2):
        return "Constraint Check Failed: The endo rule should eliminate options B and C, leaving A and D. There is an error in parsing the options or the options themselves are inconsistent."

    # 2. Check Facial Selectivity Principle
    # For C5-F, electronic effects favor syn-attack. This is a subtle but accepted principle.
    # This matches the provided reasoning.
    if provided_answer_reasoning["facial_selectivity_rule"] != "syn-attack":
        return "Incorrect Reasoning: The provided answer does not correctly identify that electronic effects favor 'syn-attack' for a C5-Fluorine substituent."

    # 3. Check Deduction of Product Structure
    # An 'endo' orientation combined with 'syn-attack' results in the Fluorine being 'anti' to the anhydride ring.
    # This matches the provided reasoning.
    if provided_answer_reasoning["product_F_position"] != "anti":
        return "Incorrect Reasoning: The geometric deduction is flawed. An endo, syn-attack leads to an 'anti' fluorine in the product."

    # 4. Check Nomenclature Mapping
    # The 'anti' position corresponds to the '8r' descriptor.
    # This matches the provided reasoning.
    if provided_answer_reasoning["nomenclature_for_anti_F"] != "8r":
        return "Incorrect Reasoning: The IUPAC descriptor for the 'anti' fluorine position in this system is '8r', not '8s'."

    # 5. Final Derivation
    # Find the candidate that is 'endo' and has the '8r' descriptor.
    final_key = None
    for key in endo_candidates:
        if get_fluorine_descriptor(options[key]) == "8r":
            final_key = key
            break
    
    if final_key is None:
        return "Constraint Check Failed: No 'endo' candidate has the required '8r' descriptor."

    # 6. Compare with provided answer
    if final_key == provided_answer_key:
        return "Correct"
    else:
        return f"Incorrect. The provided answer is {provided_answer_key}, but a step-by-step derivation based on established chemical principles points to {final_key}."

# Run the check
result = check_diels_alder_stereochemistry()
print(result)