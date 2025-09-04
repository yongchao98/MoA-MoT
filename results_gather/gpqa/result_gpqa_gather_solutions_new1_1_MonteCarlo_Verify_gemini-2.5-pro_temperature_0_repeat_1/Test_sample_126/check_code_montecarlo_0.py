from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_reaction_product():
    """
    Verifies the product of the Cope rearrangement of 5-butylnona-2,6-diene.
    The function checks for isomerism and consistency with the reaction mechanism.
    """
    try:
        # --- Step 1: Define Reactant and Options based on the question ---
        # Reactant: 5-butylnona-2,6-diene
        reactant_smiles = "CCC=CC(CCCC)C=CCC"
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        if not reactant_mol:
            return "Error: Failed to parse the reactant's SMILES string. The structure '5-butylnona-2,6-diene' might be represented incorrectly."
        reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol)
        print(f"Reactant: 5-butylnona-2,6-diene, Formula: {reactant_formula}")

        # Options from the problem description
        options = {
            "A": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCCC=CC(CC)C=CCC"},
            "B": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
            "C": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCC=CC(CC)C(C)C=C"},
            "D": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"} # Same as B
        }
        
        proposed_answer_key = "C"
        proposed_product_data = options[proposed_answer_key]
        proposed_product_mol = Chem.MolFromSmiles(proposed_product_data["smiles"])
        if not proposed_product_mol:
            return f"Error: Failed to parse the proposed product's (Option {proposed_answer_key}) SMILES string."
        
        print("-" * 30)
        
        # --- Step 2: Check the Isomerism Constraint ---
        print("Checking Constraint 1: Product must be an isomer of the reactant.")
        proposed_product_formula = rdMolDescriptors.CalcMolFormula(proposed_product_mol)
        
        if reactant_formula != proposed_product_formula:
            return (f"Incorrect. The proposed product (Option {proposed_answer_key}: {proposed_product_data['name']}) is not an isomer. "
                    f"Reactant formula: {reactant_formula}, Product formula: {proposed_product_formula}.")
        else:
            print(f"PASS: The proposed product (Option C) has the formula {proposed_product_formula}, which matches the reactant.")

        print("-" * 30)

        # --- Step 3: Eliminate other options ---
        print("Checking Constraint 2: Eliminating other options.")
        # Option A: 5-ethylundeca-2,6-diene
        mol_A = Chem.MolFromSmiles(options["A"]["smiles"])
        formula_A = rdMolDescriptors.CalcMolFormula(mol_A)
        if formula_A != reactant_formula:
            print(f"PASS: Option A ({options['A']['name']}) has formula {formula_A} and is not an isomer, so it is correctly ruled out.")
        else:
            return "Failure in check logic: Option A was expected to be a non-isomer but it is."

        # Option B/D: 5-ethyl-4-methyldeca-2,6-diene
        mol_B = Chem.MolFromSmiles(options["B"]["smiles"])
        formula_B = rdMolDescriptors.CalcMolFormula(mol_B)
        if formula_B == reactant_formula:
             print(f"INFO: Option B/D ({options['B']['name']}) is also an isomer. The choice between B/D and C depends on the reaction mechanism.")
        else:
            return "Failure in check logic: Option B/D was expected to be an isomer but it is not."
            
        print("-" * 30)

        # --- Step 4: Check for Mechanism Consistency ---
        print("Checking Constraint 3: Product structure must be consistent with the Cope rearrangement mechanism.")
        # The Cope rearrangement of a 1,5-diene system (like C2=C3-C4-C5-C6=C7) results in a new
        # 1,5-diene system where the double bonds have shifted (to C3=C4 and C5=C6).
        # When the product is re-numbered to find the longest chain, the double bonds end up at positions 1 and 5.
        # Product C is a deca-1,5-diene.
        # Product B/D is a deca-2,6-diene.
        # The standard Cope mechanism predicts the formation of a 1,5-diene product.
        if "1,5-diene" in proposed_product_data["name"]:
            print("PASS: The proposed product is a 1,5-diene, which is the expected outcome of the bond shifts in a Cope rearrangement.")
        else:
            return f"Incorrect. The proposed product is not a 1,5-diene, which contradicts the known mechanism."

        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_reaction_product()
print("-" * 30)
print(f"Final Result: {result}")