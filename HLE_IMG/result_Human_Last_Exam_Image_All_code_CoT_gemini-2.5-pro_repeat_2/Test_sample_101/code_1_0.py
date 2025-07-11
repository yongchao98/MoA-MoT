# First, ensure you have RDKit installed:
# pip install rdkit-pypi

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def solve_reaction():
    """
    Analyzes the given chemical reaction to identify Compound A.
    """
    # Define the molecules in the reaction sequence using their SMILES strings
    # Reactant 1: 3-hydroxy-pyridine-2-carbaldehyde
    smiles_reactant1 = "O=Cc1c(O)cccc1N"
    # Reactant 2: Aniline
    smiles_reactant2 = "Nc1ccccc1"
    # Intermediate: Imine formed from Step 1
    smiles_intermediate = "Oc1cnccc1C=Nc2ccccc2"
    # Reagent for Step 2: The active species are H+ and CN- (represented as HCN for mass balance)
    smiles_reagent2 = "N#C"
    # Final Product: Compound A
    smiles_product_A = "N#CC(Nc1ccccc1)c2c(O)cncc2"

    # Create RDKit molecule objects
    mol_r1 = Chem.MolFromSmiles(smiles_reactant1)
    mol_r2 = Chem.MolFromSmiles(smiles_reactant2)
    mol_intermediate = Chem.MolFromSmiles(smiles_intermediate)
    mol_reagent2 = Chem.MolFromSmiles(smiles_reagent2)
    mol_A = Chem.MolFromSmiles(smiles_product_A)

    # Calculate molecular formulas for all species
    formula_r1 = CalcMolFormula(mol_r1)
    formula_r2 = CalcMolFormula(mol_r2)
    formula_intermediate = CalcMolFormula(mol_intermediate)
    formula_reagent2 = CalcMolFormula(mol_reagent2)
    formula_A = CalcMolFormula(mol_A)
    
    # --- Output the analysis ---
    print("Reaction Analysis:")
    print("-" * 20)
    
    # Print the equation for Step 1
    print("Step 1 (Imine Formation):")
    step1_equation = f"{formula_r1} (3-hydroxy-pyridine-2-carbaldehyde) + {formula_r2} (aniline) -> {formula_intermediate} (imine) + H2O"
    print(step1_equation)

    # Print the equation for Step 2
    print("\nStep 2 (Cyanide Addition):")
    step2_equation = f"{formula_intermediate} (imine) + {formula_reagent2} (from NaCN, H+) -> {formula_A} (Compound A)"
    print(step2_equation)
    
    # Provide the identity of Compound A
    name_A = "(phenylamino)(3-hydroxypyridin-2-yl)acetonitrile"
    mw_A = Descriptors.ExactMolWt(mol_A)

    print("\n" + "-" * 20)
    print("Identity of Compound A:")
    print("-" * 20)
    print(f"Name: {name_A}")
    print(f"Molecular Formula: {formula_A}")
    print(f"Exact Molecular Weight: {mw_A:.4f} g/mol")

solve_reaction()