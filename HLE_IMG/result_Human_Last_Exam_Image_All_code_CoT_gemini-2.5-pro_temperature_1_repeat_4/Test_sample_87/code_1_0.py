# You may need to install the RDKit library first:
# pip install rdkit

import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def analyze_product(name, target_formula, proposed_smiles, explanation):
    """Analyzes a proposed product structure."""
    print(f"--- Analysis for Product {name} ---")
    print(f"Target Formula: {target_formula}")
    
    mol = Chem.MolFromSmiles(proposed_smiles)
    if mol is None:
        print("Error: Invalid SMILES string provided.")
        print("-" * 30)
        return
        
    # RDKit calculates formula without explicit hydrogens, which is what we want for comparison
    calculated_formula = CalcMolFormula(mol)
    
    print(f"Proposed Structure (SMILES): {proposed_smiles}")
    print(f"Calculated Formula: {calculated_formula}")
    
    if calculated_formula == target_formula:
        print("Formula Match: YES")
    else:
        print(f"Formula Match: NO (Calculated {calculated_formula})")
        
    print("\nProposed Origin:")
    print(explanation)
    print("-" * 40)
    print()

# Analysis for Product A
product_A_smiles = "COC(=O)[C@H]1[C@@H](OC)[C@H]2N3CCCC3N3C=CCC[C@]123"
product_A_explanation = (
    "Product A likely forms from the main cycloaddition pathway followed by further reactions:\n"
    "1. The starting material decarboxylates to form an azomethine ylide (Y).\n"
    "2. Ylide undergoes [3+2] cycloaddition with methyl propiolate to form a primary cycloadduct (P1).\n"
    "3. P1 undergoes Michael addition with methanol (CH3OH).\n"
    "4. The resulting methanol adduct is then oxidized (loses H2) to yield Product A."
)
analyze_product("A", "C14H20N2O3", product_A_smiles, product_A_explanation)

# Analysis for Product B
product_B_smiles = "CC(=O)[C@]1(C(=O)O)N(C2=NCCC2)C=CC1"
product_B_explanation = (
    "Product B likely forms from a side reaction not involving the alkyne:\n"
    "1. The starting material is deprotonated at the alpha-carbon.\n"
    "2. The resulting anion is acetylated by acetic anhydride.\n"
    "3. This intermediate (C12H16N2O3) is oxidized (loses H2 from the proline ring), creating a double bond to yield Product B."
)
analyze_product("B", "C12H14N2O3", product_B_smiles, product_B_explanation)

# Analysis for Product C
product_C_smiles = "CC(=O)[C@@]1(O)[C@H]2N(C3=NCCC3)CCC[C@@H]2C1=O"
product_C_explanation = (
    "Product C likely forms from the azomethine ylide intermediate (Y) reacting with acetic anhydride:\n"
    "1. The ylide (formed by decarboxylation) attacks acetic anhydride.\n"
    "2. This initiates a complex cascade of intramolecular rearrangement reactions.\n"
    "3. The proposed structure represents a plausible stable product from such a cascade, likely involving hydration during workup, to match the final formula."
)
analyze_product("C", "C11H16N2O3", product_C_smiles, product_C_explanation)