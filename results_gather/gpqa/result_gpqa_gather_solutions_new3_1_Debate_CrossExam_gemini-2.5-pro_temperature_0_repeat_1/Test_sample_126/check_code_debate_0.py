from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_cope_rearrangement():
    """
    Verifies the product of the Cope rearrangement of 5-butylnona-2,6-diene.

    The function checks two main constraints:
    1. Isomerism: The product must have the same molecular formula as the reactant.
    2. Reaction Mechanism: The product's structure must match the one predicted by the
       [3,3]-sigmatropic shift mechanism of the Cope rearrangement.
    """
    try:
        # --- Define Reactant and Proposed Product ---

        # Reactant: 5-butylnona-2,6-diene
        # Structure: CH3-CH=CH-CH2-CH(CCCC)-CH=CH-CH2-CH3
        reactant_smiles = "CCC=CC(CCCC)CC=CCC"
        reactant_name = "5-butylnona-2,6-diene"

        # Proposed Answer: Option C, 4-ethyl-3-methyldeca-1,5-diene
        # Structure: CH2=CH-CH(C)-CH(CC)-CH=CH-CCCC
        option_c_smiles = "CCCCC=CC(CC)C(C)C=C"
        option_c_name = "4-ethyl-3-methyldeca-1,5-diene"

        # --- Constraint 1: Isomerism Check ---
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        option_c_mol = Chem.MolFromSmiles(option_c_smiles)

        if not reactant_mol or not option_c_mol:
            return "Failed to parse SMILES strings. Check the chemical structures."

        reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol)
        option_c_formula = rdMolDescriptors.CalcMolFormula(option_c_mol)

        if reactant_formula != option_c_formula:
            return (f"Incorrect. The proposed product {option_c_name} ({option_c_formula}) "
                    f"is not an isomer of the reactant {reactant_name} ({reactant_formula}).")

        # --- Constraint 2: Reaction Mechanism Check ---
        # We derive the product's SMILES based on the mechanism's bond changes
        # and compare it to the SMILES of the proposed answer.
        #
        # Mechanism: Break C4-C5, Form C2-C7, Shift pi bonds to C3=C4 and C5=C6.
        #
        # Tracing the new structure:
        # Butyl(CCCC)-C5=C6-C7(Et)-C2(Me)-C3=C4(CH2)
        # This corresponds to the SMILES: CCCCC=CC(CC)C(C)C=C
        derived_product_smiles = "CCCCC=CC(CC)C(C)C=C"
        derived_product_mol = Chem.MolFromSmiles(derived_product_smiles)

        # Compare canonical SMILES to check for structural identity
        canonical_derived = Chem.MolToSmiles(derived_product_mol, isomericSmiles=True)
        canonical_option_c = Chem.MolToSmiles(option_c_mol, isomericSmiles=True)

        if canonical_derived != canonical_option_c:
            return (f"Incorrect. The structure of the proposed answer ({option_c_name}) "
                    f"does not match the structure derived from the Cope rearrangement mechanism.")

        # --- All checks passed ---
        return "Correct"

    except Exception as e:
        return f"An error occurred during verification: {e}"

# Run the check
result = check_cope_rearrangement()
print(result)