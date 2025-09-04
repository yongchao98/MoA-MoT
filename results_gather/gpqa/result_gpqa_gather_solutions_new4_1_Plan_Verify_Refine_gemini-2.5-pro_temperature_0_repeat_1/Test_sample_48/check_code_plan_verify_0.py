import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def parse_formula(formula_str):
    """Parses a molecular formula string (e.g., 'C8H10') into a dictionary."""
    if not formula_str:
        return {}
    pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
    matches = pattern.findall(formula_str)
    formula_dict = {}
    for element, count in matches:
        formula_dict[element] = int(count) if count else 1
    return formula_dict

def formula_to_string(formula_dict):
    """Converts a formula dictionary back to a canonical string."""
    if not formula_dict:
        return ""
    # Sort elements for canonical representation (C, H, then alphabetical)
    sorted_elements = sorted(formula_dict.keys(), key=lambda x: (x != 'C', x != 'H', x))
    
    s = ""
    for element in sorted_elements:
        count = formula_dict[element]
        s += element
        if count > 1:
            s += str(count)
    return s

def subtract_formulas(f_dict1, f_dict2):
    """Subtracts one formula dictionary from another."""
    result = f_dict1.copy()
    for element, count in f_dict2.items():
        result[element] = result.get(element, 0) - count
        if result[element] == 0:
            del result[element]
    return result

def add_formulas(f_dict1, f_dict2):
    """Adds two formula dictionaries."""
    result = f_dict1.copy()
    for element, count in f_dict2.items():
        result[element] = result.get(element, 0) + count
    return result

def get_formula_from_smiles(smiles):
    """Calculates molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return rdMolDescriptors.CalcMolFormula(mol)

def check_correctness():
    """
    Checks the correctness of the chosen answer (A) by verifying the stoichiometry
    of each reaction.
    """
    try:
        # --- Define Reactants and Products for Answer A ---
        # Reaction A: 1,1-dimethoxyethan-1-amine + but-3-en-2-ol + (H+ + Heat) ---> A
        r_a1_smiles = 'CC(N)(OC)OC'  # 1,1-dimethoxyethan-1-amine
        r_a2_smiles = 'CC(O)C=C'     # but-3-en-2-ol
        elim_a_smiles = 'CO'         # methanol (eliminated twice)
        p_a_smiles = 'NC=COC(C)=CC'  # (Z)-1-(but-2-en-2-yloxy)ethen-1-amine

        # Reaction B: (3R,4S)-3,4-dimethylhexa-1,5-diyne + Heat ---> B
        r_b_smiles = 'C#C[C@H](C)[C@@H](C)C#C' # (3R,4S)-3,4-dimethylhexa-1,5-diyne
        p_b_smiles = 'CC=C1C=CC1=CC'          # 3,4-diethylidenecyclobut-1-ene (non-stereo)

        # Reaction C: 2-((vinyloxy)methyl)but-1-ene + Heat ---> C
        r_c_smiles = 'CCC(=C)COC=C'   # 2-((vinyloxy)methyl)but-1-ene
        p_c_smiles = 'CCC(=C)CCC=O'  # 4-methylenehexanal

        # --- Verification Step 1: Reaction A Stoichiometry ---
        f_r_a1 = parse_formula(get_formula_from_smiles(r_a1_smiles))
        f_r_a2 = parse_formula(get_formula_from_smiles(r_a2_smiles))
        f_elim_a = parse_formula(get_formula_from_smiles(elim_a_smiles))
        f_p_a = parse_formula(get_formula_from_smiles(p_a_smiles))

        # Expected product formula = (Reactant1 + Reactant2) - 2 * Methanol
        expected_f_p_a = subtract_formulas(add_formulas(f_r_a1, f_r_a2), add_formulas(f_elim_a, f_elim_a))
        
        if expected_f_p_a != f_p_a:
            return (f"Incorrect: Stoichiometry for Reaction A is wrong. "
                    f"Expected product formula {formula_to_string(expected_f_p_a)}, "
                    f"but product A has formula {formula_to_string(f_p_a)}.")

        # --- Verification Step 2: Reaction B Isomerization ---
        f_r_b = get_formula_from_smiles(r_b_smiles)
        f_p_b = get_formula_from_smiles(p_b_smiles)
        if f_r_b != f_p_b:
            return (f"Incorrect: Reaction B is an isomerization, but formulas do not match. "
                    f"Reactant formula is {f_r_b}, but product B formula is {f_p_b}.")

        # --- Verification Step 3: Reaction C Isomerization ---
        f_r_c = get_formula_from_smiles(r_c_smiles)
        f_p_c = get_formula_from_smiles(p_c_smiles)
        if f_r_c != f_p_c:
            return (f"Incorrect: Reaction C is an isomerization, but formulas do not match. "
                    f"Reactant formula is {f_r_c}, but product C formula is {f_p_c}.")

        # --- Qualitative Checks (Implicitly passed if stoichiometry is correct) ---
        # C: Claisen rearrangement of allyl vinyl ether -> γ,δ-unsaturated carbonyl. 4-methylenehexanal fits.
        # B: Thermal rearrangement of 1,5-diyne -> cyclobutene derivative. Name fits.
        # A: Johnson-Claisen type reaction. Plausible intermediate.

        return "Correct"

    except (ImportError, ValueError) as e:
        return f"An error occurred during verification: {e}. Please ensure rdkit is installed."


# Run the check
result = check_correctness()
print(result)