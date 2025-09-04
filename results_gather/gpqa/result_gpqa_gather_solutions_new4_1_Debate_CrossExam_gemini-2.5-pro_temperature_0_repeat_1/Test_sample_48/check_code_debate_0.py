from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_formula_from_smiles(smiles: str) -> str:
    """Calculates the molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return rdMolDescriptors.CalcMolFormula(mol)
    return None

def check_correctness():
    """
    Checks the correctness of the final answer for the given organic chemistry problem.
    The function verifies the products of option 'C' based on stoichiometry and
    known reaction mechanisms.
    """
    # The final answer provided by the LLM aggregator to be checked.
    final_answer = 'C'

    # A database of SMILES strings for all relevant reactants and products.
    smiles_db = {
        # Reactants
        "1,1-dimethoxyethan-1-amine": "CC(N)(OC)OC",
        "but-3-en-2-ol": "C=CC(C)O",
        "(3R,4S)-3,4-dimethylhexa-1,5-diyne": "C#C[C@H](C)[C@H](C)C#C",
        "2-((vinyloxy)methyl)but-1-ene": "C=C(CC)CO-C=C",

        # Products for A
        "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine": "C/C=C(O/C(N)=C)/C",
        "6-methyl-3,4-dihydro-2H-pyran-2-amine": "CC1=CCCC(N)O1",

        # Products for B
        "(3Z,4E)-3,4-diethylidenecyclobut-1-ene": "C/C=C1/C(=C/C)C=C1",
        "(1Z,2E)-1,2-diethylidenecyclobutane": "C/C=C1/C(=C/C)CCC1",

        # Products for C
        "4-methylenehexanal": "O=CCCC(=C)CC",
        "4-methylenehexan-1-ol": "OCCCC(=C)CC",
    }

    # The product names for each option as given in the question.
    products_map = {
        "A": {
            "A_name": "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine",
            "B_name": "(3Z,4E)-3,4-diethylidenecyclobut-1-ene",
            "C_name": "4-methylenehexan-1-ol"
        },
        "B": {
            "A_name": "6-methyl-3,4-dihydro-2H-pyran-2-amine",
            "B_name": "(1Z,2E)-1,2-diethylidenecyclobutane",
            "C_name": "4-methylenehexanal"
        },
        "C": {
            "A_name": "(Z)-1-(but-2-en-2-yloxy)ethen-1-amine",
            "B_name": "(3Z,4E)-3,4-diethylidenecyclobut-1-ene",
            "C_name": "4-methylenehexanal"
        },
        "D": {
            "A_name": "6-methyl-3,4-dihydro-2H-pyran-2-amine",
            "B_name": "(1Z,2E)-1,2-diethylidenecyclobutane",
            "C_name": "4-methylenehexan-1-ol"
        }
    }

    # --- Check Reaction 1 (A) ---
    # This is a condensation reaction eliminating two molecules of methanol (2 * CH4O).
    # Reactants: C4H11NO2 + C4H8O = C8H19NO3
    # Product: C8H19NO3 - C2H8O2 = C6H11NO
    expected_p1_formula = "C6H11NO"
    
    product_A_name = products_map[final_answer]["A_name"]
    product_A_formula = get_formula_from_smiles(smiles_db[product_A_name])

    if product_A_formula != expected_p1_formula:
        return (f"Incorrect: Product A ('{product_A_name}') has formula {product_A_formula}, "
                f"but the expected formula from the condensation reaction is {expected_p1_formula}.")

    # --- Check Reaction 2 (B) ---
    # This is a thermal rearrangement (isomerization). Product formula must match reactant formula.
    r2_formula = get_formula_from_smiles(smiles_db["(3R,4S)-3,4-dimethylhexa-1,5-diyne"]) # C8H10
    
    product_B_name = products_map[final_answer]["B_name"]
    product_B_formula = get_formula_from_smiles(smiles_db[product_B_name])

    if product_B_formula != r2_formula:
        return (f"Incorrect: Product B ('{product_B_name}') has formula {product_B_formula}, "
                f"but it should be an isomer of the reactant with formula {r2_formula}.")
    
    # Mechanistic check: The Hopf rearrangement of a 1,5-diyne yields a cyclobutene, not a cyclobutane.
    if "cyclobutane" in product_B_name:
        return (f"Incorrect: Product B ('{product_B_name}') is a cyclobutane. The thermal rearrangement "
                f"of a 1,5-diyne should yield a cyclobutene derivative.")

    # --- Check Reaction 3 (C) ---
    # This is a Claisen rearrangement (isomerization). Product formula must match reactant formula.
    r3_formula = get_formula_from_smiles(smiles_db["2-((vinyloxy)methyl)but-1-ene"]) # C7H12O
    
    product_C_name = products_map[final_answer]["C_name"]
    product_C_formula = get_formula_from_smiles(smiles_db[product_C_name])

    if product_C_formula != r3_formula:
        return (f"Incorrect: Product C ('{product_C_name}') has formula {product_C_formula}, "
                f"but it should be an isomer of the reactant with formula {r3_formula}.")

    # Mechanistic check: Claisen rearrangement of an allyl vinyl ether yields a carbonyl, not an alcohol.
    if "ol" in product_C_name:
        return (f"Incorrect: Product C ('{product_C_name}') is an alcohol. The Claisen rearrangement "
                f"should yield a carbonyl compound (aldehyde or ketone).")

    # If all checks pass for the given answer 'C'.
    return "Correct"

# Execute the check and print the result.
result = check_correctness()
print(result)