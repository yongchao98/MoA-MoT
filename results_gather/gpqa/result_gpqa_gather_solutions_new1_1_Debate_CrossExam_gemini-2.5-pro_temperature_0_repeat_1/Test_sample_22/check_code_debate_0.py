from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def get_molecular_formula(smiles: str) -> str:
    """Calculates the molecular formula from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES"
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        return f"Error: {e}"

def check_answer():
    """
    Checks the correctness of the answer to the organic chemistry question.
    """
    # --- Define Reactants and Products using SMILES strings ---
    reactant_smiles = "c1ccc(cc1)OCC(C)(C)C=C"  # ((2,2-dimethylbut-3-en-1-yl)oxy)benzene

    # Option A: 2-(2,2-dimethylbutyl)phenol and 4-(2,2-dimethylbutyl)phenol
    option_a_smiles = ["c1ccc(c(c1)O)C(C)(C)CC", "c1cc(ccc1O)C(C)(C)CC"]

    # Option B: 3,3,4-trimethylchromane and 3-isopropyl-3-methyl-2,3-dihydrobenzofuran
    option_b_smiles = ["CC1(C)C(C)C2=CC=CC=C2OC1", "CC(C)C1(C)OC2=CC=CC=C2C1"]

    # Option C: (4-bromo-2,2-dimethylbutoxy)benzene and (3-bromo-2,2-dimethylbutoxy)benzene
    option_c_smiles = ["c1ccc(cc1)OCC(C)(C)CCBr", "c1ccc(cc1)OCC(C)(C)C(C)Br"]
    
    # Option D: (4-bromo-2,2-dimethylbutoxy)benzene and ((2,3-dimethylbut-2-en-1-yl)oxy)benzene
    option_d_smiles = ["c1ccc(cc1)OCC(C)(C)CCBr", "c1ccc(cc1)OCC(=C(C)C)C"]

    # The final answer provided by the analysis
    final_answer = "B"

    # --- Perform Checks ---
    reactant_formula = get_molecular_formula(reactant_smiles)
    
    # Check Option A
    formulas_a = [get_molecular_formula(s) for s in option_a_smiles]
    # Expected transformation: Ether cleavage and reduction (addition of H2).
    # HBr is an acid/nucleophile source, not a reducing agent.
    if any(f != "C12H18O" for f in formulas_a):
         return f"Incorrect: Option A has incorrect molecular formulas. Expected C12H18O, but got {formulas_a}."
    if "C12H18O" in formulas_a:
        return "Incorrect: The products in Option A (C12H18O) would require reduction (saturation of the double bond), which is not a reaction caused by HBr. The reactant is C12H16O."

    # Check Option C
    formulas_c = [get_molecular_formula(s) for s in option_c_smiles]
    # Expected transformation: Addition of HBr. Formula should be C12H16O + HBr -> C12H17BrO
    if any(f != "C12H17BrO" for f in formulas_c):
        return f"Incorrect: Option C is an addition reaction, so products should have the formula C12H17BrO. Instead, got {formulas_c}."
    # This check passes, but the mechanism is less likely than cyclization for explaining two major products.

    # Check Option D
    formulas_d = [get_molecular_formula(s) for s in option_d_smiles]
    # This option mixes two different reaction types.
    if formulas_d[0] != "C12H17BrO" or formulas_d[1] != "C12H16O":
        return f"Incorrect: Option D has inconsistent molecular formulas. Expected one addition product (C12H17BrO) and one isomer (C12H16O), but got {formulas_d}."
    if formulas_d[0] == "C12H17BrO" and formulas_d[1] == "C12H16O":
        return "Incorrect: Option D proposes that the two main products come from two completely different reaction types (addition and isomerization). This is mechanistically inconsistent for explaining the two observed spots."

    # Check Option B
    formulas_b = [get_molecular_formula(s) for s in option_b_smiles]
    # Expected transformation: Intramolecular cyclization, which is an isomerization.
    # Product formulas must match the reactant formula.
    if any(f != reactant_formula for f in formulas_b):
        return f"Incorrect: Option B should represent an isomerization, so product formulas must match the reactant ({reactant_formula}). Instead, got {formulas_b}."
    
    # If all checks for B pass, it's the most plausible answer.
    if final_answer == "B":
        return "Correct"
    else:
        return f"Incorrect: The provided answer is {final_answer}, but the analysis shows B is the only chemically plausible option. The reaction is an intramolecular cyclization, which is an isomerization. Only the products in Option B are isomers of the reactant ({reactant_formula})."

# Run the check
result = check_answer()
print(result)