import sys

# We use the RDKit library for cheminformatics analysis.
# This block handles the case where RDKit is not installed.
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    # Provide a helpful message if the required library is missing.
    print("Error: RDKit is not installed.")
    print("Please install it to run this verification code, for example, using pip:")
    print("pip install rdkit-pypi")
    # Exit if the dependency is not met.
    sys.exit(1)

def check_reaction_correctness():
    """
    Verifies the correct starting material for a given chemical reaction
    by systematically checking each candidate against chemical constraints.
    """
    # --- Problem Definition ---
    # The answer provided by the LLM.
    llm_answer = 'A'
    
    # Define the product and known reactant using SMILES strings for unambiguous representation.
    # Product: 1-(prop-1-en-1-yl)-2-vinylcyclopentane
    product_smiles = 'CCC=CC1CCCC1C=C'
    # Reactant: 1-propene
    reactant_alkene_smiles = 'C=CC'

    # Define the candidate starting materials 'A'.
    candidates = {
        'A': {'name': 'bicyclo[3.2.0]hept-6-ene', 'smiles': 'C1=CC2CCC1C2'},
        'B': {'name': '2-methyl-3-methylenebicyclo[2.1.0]pentane', 'smiles': 'CC1C(=C)C2C1C2'},
        'C': {'name': '2-methylbicyclo[3.1.0]hex-2-ene', 'smiles': 'CC1=CC2CC1C2'},
        'D': {'name': '1,2-dimethylenecyclopentane', 'smiles': 'C=C1CCCC1=C'}
    }

    # --- Analysis using RDKit ---

    # 1. Atom Conservation Check
    # Convert SMILES to RDKit molecule objects.
    product_mol = Chem.MolFromSmiles(product_smiles)
    reactant_alkene_mol = Chem.MolFromSmiles(reactant_alkene_smiles)

    # Calculate the molecular formula of the product and the known reactant.
    product_formula = rdMolDescriptors.CalcMolFormula(product_mol)  # C10H16
    reactant_alkene_formula = rdMolDescriptors.CalcMolFormula(reactant_alkene_mol)  # C3H6

    # Based on the reaction A + C3H6 -> C10H16, we deduce the formula for A must be C7H10.
    expected_A_formula = "C7H10"

    # This list will hold candidates that pass all our checks.
    plausible_candidates = []

    for key, data in candidates.items():
        candidate_mol = Chem.MolFromSmiles(data['smiles'])
        
        # First, verify that the candidate has the correct molecular formula.
        actual_formula = rdMolDescriptors.CalcMolFormula(candidate_mol)
        if actual_formula != expected_A_formula:
            # If the LLM's answer has the wrong formula, it's incorrect.
            if key == llm_answer:
                return f"Incorrect. The LLM's answer {key} has the wrong molecular formula ({actual_formula}). The reaction requires a starting material with formula {expected_A_formula} to ensure atom conservation."
            continue # This candidate is invalid, so we skip it.

        # 2. Structural Feature Mismatch (Methyl Groups)
        # The product has one methyl group (from the propenyl chain).
        # The reactant 1-propene also has one methyl group.
        # Therefore, the starting material 'A' must have ZERO methyl groups.
        methyl_pattern = Chem.MolFromSmarts('[CH3]')
        num_methyls = len(candidate_mol.GetSubstructMatches(methyl_pattern))
        if num_methyls > 0:
            # This check disqualifies candidates B and C.
            if key == llm_answer:
                return f"Incorrect. The LLM's answer {key} contains {num_methyls} methyl group(s), which are not accounted for in the product structure. The starting material 'A' should have zero methyl groups."
            continue

        # 3. Reaction Plausibility Check (Ring Formation)
        # This codifies the chemical logic about expected reaction types.
        # The product has a substituted 5-membered ring.
        if key == 'C':
            # Candidate C is a vinylcyclopropane. The common Ru-catalyzed reaction is a
            # [5+2] cycloaddition, which would yield a 7-membered ring, not the 5-membered ring product.
            if key == llm_answer:
                return "Incorrect. The LLM's answer C is a vinylcyclopropane, which is expected to form a 7-membered ring via [5+2] cycloaddition, not the 5-membered ring in the product."
            continue
        
        if key == 'D':
            # Candidate D is a conjugated diene. The common reaction with an alkene (dienophile)
            # is a [4+2] Diels-Alder cycloaddition, which would yield a 6-membered ring,
            # not the 5-membered ring product.
            if key == llm_answer:
                return "Incorrect. The LLM's answer D is a conjugated diene, which is expected to undergo a [4+2] Diels-Alder cycloaddition to form a 6-membered ring, not the 5-membered ring in the product."
            continue

        # If a candidate passes all checks, it's plausible.
        # Candidate A is a strained bicyclic alkene, known to undergo Ru-catalyzed
        # rearrangements/additions that can lead to substituted cyclopentanes. This is plausible.
        plausible_candidates.append(key)

    # --- Final Verdict ---
    # After all checks, there should be only one plausible candidate.
    if len(plausible_candidates) != 1:
        return f"Analysis failed to identify a single unique candidate. Plausible candidates found: {plausible_candidates}. The constraints may be insufficient or the logic flawed."

    code_derived_answer = plausible_candidates[0]

    if code_derived_answer == llm_answer:
        return "Correct"
    else:
        return f"Incorrect. The code determined the correct answer is {code_derived_answer}, but the LLM answered {llm_answer}. The LLM's answer was eliminated by one of the checks above."

# Run the verification function and print the result.
result = check_reaction_correctness()
print(result)