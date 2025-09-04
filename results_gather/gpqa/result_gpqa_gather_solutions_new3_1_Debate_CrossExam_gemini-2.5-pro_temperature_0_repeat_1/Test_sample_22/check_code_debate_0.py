from rdkit import Chem
from rdkit.Chem import Descriptors

def check_chemistry_answer():
    """
    Checks the correctness of the answer by verifying the chemical principles.
    The core principle here is that an intramolecular cyclization is an isomerization reaction,
    so the products must have the same molecular formula as the reactant.
    """
    # The final answer provided by the LLM analysis
    final_answer = 'A'

    # Define the reactant and products from the options using SMILES strings
    # SMILES (Simplified Molecular Input Line Entry System) is a standard way to represent molecules.
    reactant_smiles = 'c1ccccc1OCC(C)(C)C=C'  # ((2,2-dimethylbut-3-en-1-yl)oxy)benzene

    options = {
        'A': {
            "names": ["3,3,4-trimethylchromane", "3-isopropyl-3-methyl-2,3-dihydrobenzofuran"],
            "smiles": [
                'c1ccc2c(c1)OC(C(C)C)C2(C)C',
                'c1ccc2c(c1)OC(C(C)C)(C)C2'
            ]
        },
        'B': {
            "names": ["2-(2,2-dimethylbutyl)phenol", "4-(2,2-dimethylbutyl)phenol"],
            "smiles": [
                'c1c(cccc1O)CC(C)(C)CC',
                'c1cc(ccc1O)CC(C)(C)CC'
            ]
        },
        'C': {
            "names": ["(4-bromo-2,2-dimethylbutoxy)benzene", "(3-bromo-2,2-dimethylbutoxy)benzene"],
            "smiles": [
                'c1ccccc1OCC(C)(C)CCBr',
                'c1ccccc1OCC(C)(C)C(C)Br'
            ]
        },
        'D': {
            "names": ["(4-bromo-2,2-dimethylbutoxy)benzene", "((2,3-dimethylbut-2-en-1-yl)oxy)benzene"],
            "smiles": [
                'c1ccccc1OCC(C)(C)CCBr',
                'c1ccccc1OCC(=C(C)C)C'
            ]
        }
    }

    # --- Verification Logic ---
    try:
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        if not reactant_mol:
            return "Error: Could not parse the reactant SMILES string."
        reactant_formula = Descriptors.rdMolDescriptors.CalcMolFormula(reactant_mol)
    except Exception as e:
        return f"An error occurred with rdkit: {e}. Please ensure rdkit is installed (`pip install rdkit`)."

    # Check the provided final answer
    selected_option_data = options.get(final_answer)
    if not selected_option_data:
        return f"Error: The final answer '{final_answer}' is not a valid option."

    # Check if both products in the selected option are isomers of the reactant
    product1_mol = Chem.MolFromSmiles(selected_option_data["smiles"][0])
    product2_mol = Chem.MolFromSmiles(selected_option_data["smiles"][1])
    
    if not product1_mol or not product2_mol:
        return f"Error: Could not parse product SMILES strings for option {final_answer}."

    product1_formula = Descriptors.rdMolDescriptors.CalcMolFormula(product1_mol)
    product2_formula = Descriptors.rdMolDescriptors.CalcMolFormula(product2_mol)

    if product1_formula == reactant_formula and product2_formula == reactant_formula:
        # The answer is consistent with the principle of isomerization.
        # Now, let's confirm why other options are wrong.
        error_messages = []
        for key, data in options.items():
            if key == final_answer:
                continue
            p1_formula = Descriptors.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(data["smiles"][0]))
            p2_formula = Descriptors.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(data["smiles"][1]))
            if p1_formula != reactant_formula or p2_formula != reactant_formula:
                # This confirms other options are not simple isomerizations.
                pass
            else:
                # This case should not be hit if A is uniquely correct.
                error_messages.append(f"Logic error: Option {key} also contains isomers.")
        
        if not error_messages:
            return "Correct"
        else:
            return "\n".join(error_messages)
    else:
        # The selected answer is incorrect because its products are not isomers.
        reason = f"Incorrect. The products in option {final_answer} are not isomers of the reactant.\n"
        reason += f"Reactant formula: {reactant_formula}\n"
        reason += f"Product 1 ('{selected_option_data['names'][0]}') formula: {product1_formula}\n"
        reason += f"Product 2 ('{selected_option_data['names'][1]}') formula: {product2_formula}\n"
        
        if "Br" in product1_formula or "Br" in product2_formula:
            reason += "This indicates an addition of HBr, not an intramolecular cyclization."
        elif "H" in product1_formula and int(product1_formula.split('H')[1].split('O')[0]) > int(reactant_formula.split('H')[1].split('O')[0]):
             reason += "This indicates a reduction reaction (addition of H2), which is not possible with HBr."
        else:
             reason += "The reaction should be an isomerization."
        return reason

# Execute the check
print(check_chemistry_answer())