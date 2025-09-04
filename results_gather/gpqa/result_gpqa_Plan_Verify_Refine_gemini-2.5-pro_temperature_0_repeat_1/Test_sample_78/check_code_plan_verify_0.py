import sys
from io import StringIO

def check_answer():
    """
    Checks the correctness of the provided LLM answer by:
    1. Verifying the deduced product structure against the problem's constraints.
    2. Verifying the chosen reactant (Compound X) is a valid precursor.
    3. Verifying the reasons for eliminating other options.
    """
    # Suppress RDKit's non-critical error messages for cleaner output
    original_stderr = sys.stderr
    sys.stderr = StringIO()
    try:
        from rdkit import Chem
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    except ImportError:
        return "Execution Error: RDKit is not installed. Please install it using 'pip install rdkit-pypi' to run this verification code."
    finally:
        # Restore stderr
        sys.stderr = original_stderr

    # --- Definitions ---
    # Molecular formula given in the question for Compound X
    target_formula = "C11H12O"
    
    # The final answer provided by the LLM
    llm_answer = 'A'

    # SMILES representations of the compounds involved
    structures = {
        'Product': {
            'name': '4-(p-tolyl)but-3-en-2-one',
            'smiles': 'CC(=O)C=Cc1ccc(C)cc1'
        },
        'A': {
            'name': '2-(4-methylstyryl)oxirane',
            'smiles': 'Cc1ccc(C=CC2OC2)cc1'
        },
        'B': {
            'name': '2-methyl-3-styryloxirane',
            'smiles': 'CC1OC1C=Cc1ccccc1'
        },
        'C': {
            'name': '2-styrylepoxide',
            'smiles': 'c1ccccc1C=CC1OC1'
        },
        'D': {
            'name': '2-(1-phenylprop-1-en-2-yl)oxirane',
            'smiles': 'C(=C(c1ccccc1)C)C1OC1'
        }
    }

    # Helper functions
    def get_formula(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return CalcMolFormula(mol) if mol else "Invalid SMILES"

    def has_substructure(smiles, smarts):
        mol = Chem.MolFromSmiles(smiles)
        pattern = Chem.MolFromSmarts(smarts)
        return mol.HasSubstructMatch(pattern) if mol and pattern else False

    # --- Verification Logic ---
    errors = []

    # 1. Verify the deduced product structure
    product_info = structures['Product']
    product_formula = get_formula(product_info['smiles'])
    if product_formula != target_formula:
        errors.append(f"Reasoning Error: The deduced product '{product_info['name']}' has formula {product_formula}, but the reactant's formula is {target_formula}.")
    
    # Check for ketone group (matches 13C NMR δ 197.7)
    if not has_substructure(product_info['smiles'], '[#6]C(=O)[#6]'):
        errors.append(f"Reasoning Error: The deduced product structure lacks a ketone, which contradicts the 13C NMR signal at 197.7 ppm.")
        
    # Check for p-tolyl group (matches 1H NMR para-substitution pattern)
    if not has_substructure(product_info['smiles'], 'c1(C)ccc(cc1)'):
        errors.append(f"Reasoning Error: The deduced product structure lacks a p-tolyl group, which contradicts the 1H NMR A2B2 pattern.")

    # 2. Verify the chosen answer (Option A)
    reactant_info = structures[llm_answer]
    reactant_formula = get_formula(reactant_info['smiles'])
    
    if reactant_formula != target_formula:
        errors.append(f"Constraint Violation: The chosen answer, Option {llm_answer} ('{reactant_info['name']}'), has formula {reactant_formula}, which does not match the required {target_formula}.")
    
    if not has_substructure(reactant_info['smiles'], 'c1(C)ccc(cc1)'):
        errors.append(f"Plausibility Error: The chosen reactant (Option {llm_answer}) must contain a p-tolyl group to form the product, but it does not.")

    # 3. Verify the elimination of other options based on the LLM's reasoning
    # Option B: Phenyl group, not p-tolyl.
    option_b_info = structures['B']
    if has_substructure(option_b_info['smiles'], 'c1(C)ccc(cc1)'):
        errors.append(f"Analysis Error: The LLM correctly states Option B has a phenyl group, but the structure '{option_b_info['smiles']}' contains a p-tolyl group.")
    if not has_substructure(option_b_info['smiles'], 'c1ccccc1'):
        errors.append(f"Analysis Error: The LLM states Option B has a phenyl group, but the structure '{option_b_info['smiles']}' does not.")

    # Option C: Wrong formula (C10H10O)
    option_c_info = structures['C']
    option_c_formula = get_formula(option_c_info['smiles'])
    if option_c_formula == target_formula:
        errors.append(f"Analysis Error: The LLM claims Option C has the wrong formula, but its formula ({option_c_formula}) matches the target.")
    if option_c_formula != 'C10H10O':
        errors.append(f"Analysis Error: The LLM's reason for dismissing Option C (formula is C10H10O) is incorrect; its formula is actually {option_c_formula}.")

    # Option D: Phenyl group, not p-tolyl.
    option_d_info = structures['D']
    if has_substructure(option_d_info['smiles'], 'c1(C)ccc(cc1)'):
        errors.append(f"Analysis Error: The LLM correctly states Option D has a phenyl group, but the structure '{option_d_info['smiles']}' contains a p-tolyl group.")
    if not has_substructure(option_d_info['smiles'], 'c1ccccc1'):
        errors.append(f"Analysis Error: The LLM states Option D has a phenyl group, but the structure '{option_d_info['smiles']}' does not.")

    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. The following issues were found:\n- " + "\n- ".join(errors)

# Run the check and print the result
result = check_answer()
print(result)