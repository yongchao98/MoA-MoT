from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_cope_rearrangement_answer():
    """
    Checks the correctness of the answer for the Cope rearrangement of 5-butylnona-2,6-diene.
    """
    # --- Step 1: Define Reactant and Options ---
    reactant_name = "5-butylnona-2,6-diene"
    # SMILES for 5-butylnona-2,6-diene: CH3-CH=CH-CH2-CH(C4H9)-CH=CH-CH2-CH3
    reactant_smiles = "CCC=CC(CCCC)CC=CC"
    
    options = {
        "A": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCC=CC(CC)C(C)C=C"},
        "B": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
        "C": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCC=CC(CC)CC=CC"},
        "D": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"} # Identical to B
    }
    
    llm_answer = "B"

    # --- Step 2: Verify Isomerism Constraint ---
    try:
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol)
    except:
        return "Error: Could not parse the reactant SMILES string."

    for key, data in options.items():
        if key == "D": continue # Skip duplicate
        try:
            mol = Chem.MolFromSmiles(data["smiles"])
            formula = rdMolDescriptors.CalcMolFormula(mol)
            if formula != reactant_formula:
                return f"Incorrect. The reaction is an isomerization, but option {key} ({data['name']}) has formula {formula}, which does not match the reactant's formula {reactant_formula}."
        except:
            return f"Error: Could not parse the SMILES string for option {key}."

    # --- Step 3: Determine the Correct Product via Chemical Principles ---
    # The Cope rearrangement of 5-butylnona-2,6-diene is a well-defined [3,3]-sigmatropic shift.
    # The 1,5-diene system involves carbons C2 through C7 of the nonane chain.
    # 1. The C4-C5 sigma bond breaks.
    # 2. A new C2-C7 sigma bond forms.
    # 3. The pi bonds shift to C3=C4 and C5=C6.
    # Tracing the atoms through this rearrangement results in a new molecule.
    # The IUPAC name of this resulting molecule is 4-ethyl-3-methyldeca-1,5-diene.
    
    correct_answer_key = "A"
    
    # --- Step 4: Compare with LLM's Answer and Return Verdict ---
    if llm_answer == correct_answer_key:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer is {llm_answer}, but the correct product is A.\n"
            f"The reaction is a Cope rearrangement, a [3,3]-sigmatropic shift that occurs when a 1,5-diene is heated.\n"
            f"1. Reactant: {reactant_name} is a 1,5-diene system.\n"
            f"2. Mechanism: The bond between C4 and C5 breaks, a new bond between C2 and C7 forms, and the double bonds shift accordingly.\n"
            f"3. Product: This rearrangement transforms the carbon skeleton into a new structure. The correct IUPAC name for this product is {options[correct_answer_key]['name']}.\n"
            f"4. Conclusion: The correct option is A. Option {llm_answer} ({options[llm_answer]['name']}) is an isomer but is not the product of this specific rearrangement."
        )
        return reason

# Execute the check
result = check_cope_rearrangement_answer()
print(result)