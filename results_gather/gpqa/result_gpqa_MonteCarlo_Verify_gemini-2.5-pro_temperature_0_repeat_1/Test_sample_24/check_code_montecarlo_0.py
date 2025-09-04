import sys
from collections import defaultdict
import re

# This script requires the RDKit library.
# You can install it via pip: pip install rdkit
# Or via conda: conda install -c conda-forge rdkit

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("Error: RDKit library not found.")
    print("Please install it using 'pip install rdkit' or 'conda install -c conda-forge rdkit'")
    sys.exit(1)

def parse_formula(formula: str) -> defaultdict:
    """Parses a molecular formula string into a dictionary of element counts."""
    counts = defaultdict(int)
    for element, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        counts[element] += int(count) if count else 1
    return counts

def check_reaction_1(reactant_name: str, product_name: str) -> (bool, str):
    """
    Checks the validity of the proposed reactants for Reaction 1.
    Reaction: A + H2SO4 ---> 2,8-dimethylspiro[4.5]decan-6-one
    This is a Pinacol rearrangement, which converts a 1,2-diol to a ketone with acid.
    """
    # --- Define Structures ---
    # Reactant A from option C: 2,7-dimethyloctahydronaphthalene-4a,8a-diol
    # This is a decalin system with a 1,2-diol at the bridgehead carbons.
    # A representative SMILES for a valid stereoisomer that undergoes the rearrangement:
    reactant_smiles = "CC1CC[C@]2(O)CCCC(C)C[C@]12O"
    
    # Product: 2,8-dimethylspiro[4.5]decan-6-one
    product_smiles = "CC1CCC(C)C2(C1)CCCC2=O"

    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    if not reactant_mol or not product_mol:
        return False, f"Internal error: Failed to parse SMILES for Reaction 1. Reactant: {reactant_smiles}, Product: {product_smiles}"

    # --- Check Constraints ---
    # Constraint 1: Reactant must be a 1,2-diol (vicinal diol).
    diol_pattern = Chem.MolFromSmarts("[#6]([#8H1])-[#6]([#8H1])")
    if not reactant_mol.HasSubstructMatch(diol_pattern):
        return False, f"Constraint failed for Reaction 1: Reactant '{reactant_name}' is not a 1,2-diol, which is required for a Pinacol rearrangement."

    # Constraint 2: Product must be a ketone.
    ketone_pattern = Chem.MolFromSmarts("[#6][C](=O)[#6]")
    if not product_mol.HasSubstructMatch(ketone_pattern):
        return False, f"Constraint failed for Reaction 1: Product '{product_name}' is not a ketone."

    # Constraint 3: The reaction is a dehydration. Formula of reactant should be formula of product + H2O.
    reactant_formula = Descriptors.rdMolDescriptors.CalcMolFormula(reactant_mol)
    product_formula = Descriptors.rdMolDescriptors.CalcMolFormula(product_mol)
    
    r_counts = parse_formula(reactant_formula)
    p_counts = parse_formula(product_formula)
    
    # Expected reactant counts based on product + H2O
    expected_r_counts = p_counts.copy()
    expected_r_counts['H'] += 2
    expected_r_counts['O'] += 1

    if r_counts != expected_r_counts:
        # Create expected formula string for a more readable error message
        expected_formula_str = "".join([f"{el}{cnt}" for el, cnt in sorted(expected_r_counts.items())])
        return False, f"Constraint failed for Reaction 1: Stoichiometry is incorrect. Reactant formula is {reactant_formula}, but expected {expected_formula_str} (Product {product_formula} + H2O)."

    return True, "Reaction 1 is a valid Pinacol rearrangement."

def check_reaction_2(reactant_name: str, product_name: str) -> (bool, str):
    """
    Checks the validity of the proposed reactants for Reaction 2.
    Reaction: B + BuLi + H+ ---> 4-methyl-1-phenylpent-3-en-1-ol
    This is a Wittig rearrangement, which is an isomerization of an ether to an alcohol.
    """
    # --- Define Structures ---
    # Reactant B from option C: (((3-methylbut-2-en-1-yl)oxy)methyl)benzene (benzyl prenyl ether)
    reactant_smiles = "CC(C)=CCOCc1ccccc1"
    
    # Product: 4-methyl-1-phenylpent-3-en-1-ol
    # Structure: Ph-CH(OH)-CH2-CH=C(CH3)2
    product_smiles = "CC(C)=CHCH2C(O)c1ccccc1"

    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    if not reactant_mol or not product_mol:
        return False, f"Internal error: Failed to parse SMILES for Reaction 2. Reactant: {reactant_smiles}, Product: {product_smiles}"

    # --- Check Constraints ---
    # Constraint 1: Reactant must be an ether.
    ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")
    if not reactant_mol.HasSubstructMatch(ether_pattern):
        return False, f"Constraint failed for Reaction 2: Reactant '{reactant_name}' is not an ether."

    # Constraint 2: Product must be an alcohol.
    alcohol_pattern = Chem.MolFromSmarts("[#6]-[OH1]")
    if not product_mol.HasSubstructMatch(alcohol_pattern):
        return False, f"Constraint failed for Reaction 2: Product '{product_name}' is not an alcohol."

    # Constraint 3: The reaction is a rearrangement (isomerization). Molecular formulas must be identical.
    reactant_formula = Descriptors.rdMolDescriptors.CalcMolFormula(reactant_mol)
    product_formula = Descriptors.rdMolDescriptors.CalcMolFormula(product_mol)

    if reactant_formula != product_formula:
        return False, f"Constraint failed for Reaction 2: It should be an isomerization, but formulas differ. Reactant: {reactant_formula}, Product: {product_formula}."

    # Constraint 4: The specific rearrangement must be chemically correct.
    # Deprotonation at the benzylic position leads to a [1,2]-Wittig rearrangement.
    # Reactant: (CH3)2C=CH-CH2-O-CH2-Ph -> Product: (CH3)2C=CH-CH2-CH(OH)-Ph
    # The structure of this predicted product matches the given product.
    # We can verify this by comparing canonical SMILES.
    
    pred_product_mol = Chem.MolFromSmiles("CC(C)=CHCH2C(O)c1ccccc1")
    
    pred_canon_smiles = Chem.MolToSmiles(pred_product_mol, canonical=True)
    given_canon_smiles = Chem.MolToSmiles(product_mol, canonical=True)

    if pred_canon_smiles != given_canon_smiles:
        return False, f"Constraint failed for Reaction 2: The rearrangement of reactant B does not yield the specified product."

    return True, "Reaction 2 is a valid [1,2]-Wittig rearrangement."

def check_answer():
    """
    Main function to check the correctness of the LLM's answer.
    """
    # The LLM's final answer is 'C'.
    llm_answer_key = 'C'
    
    options = {
        'C': {
            'A': '2,7-dimethyloctahydronaphthalene-4a,8a-diol',
            'B': '(((3-methylbut-2-en-1-yl)oxy)methyl)benzene'
        }
    }
    
    # The question's products
    product_1_name = "2,8-dimethylspiro[4.5]decan-6-one"
    product_2_name = "4-methyl-1-phenylpent-3-en-1-ol"

    selected_option = options.get(llm_answer_key)
    if not selected_option:
        print(f"The provided answer '{llm_answer_key}' is not one of the options.")
        return

    reactant_A_name = selected_option['A']
    reactant_B_name = selected_option['B']

    # Check Reaction 1
    is_correct_1, reason_1 = check_reaction_1(reactant_A_name, product_1_name)
    if not is_correct_1:
        print(f"The provided answer '{llm_answer_key}' is incorrect.")
        print(f"Reason: {reason_1}")
        return

    # Check Reaction 2
    is_correct_2, reason_2 = check_reaction_2(reactant_B_name, product_2_name)
    if not is_correct_2:
        print(f"The provided answer '{llm_answer_key}' is incorrect.")
        print(f"Reason: {reason_2}")
        return

    # If both checks pass, the answer is correct.
    print("Correct")

# Execute the check
check_answer()