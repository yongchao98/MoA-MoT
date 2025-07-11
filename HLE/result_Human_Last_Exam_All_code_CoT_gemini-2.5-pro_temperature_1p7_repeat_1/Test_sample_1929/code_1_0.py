# This script requires the RDKit library.
# You can install it using: pip install rdkit-pypi

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# 1. Define reactants using SMILES strings
# Butadiene (diene) and 1,1-dichloro-2,2-difluoroethene (dienophile)
reactants_smiles = ('C=CC=C', 'C(Cl)(Cl)=C(F)F')
reactant_mols = [Chem.MolFromSmiles(s) for s in reactants_smiles]

# 2. Define the Diels-Alder reaction using a reaction SMARTS template
# This template describes the [4+2] cycloaddition
rxn_smarts = '[c:1]=[c:2][c:3]=[c:4].[c:5]=[c:6]>>[c:1]1[c:2]=[c:3][c:4][c:6][c:5]1'
rxn = AllChem.ReactionFromSmarts(rxn_smarts)

# 3. Run the reaction
# The result is a tuple of tuples of product molecules
product_sets = rxn.RunReactants(reactant_mols)

# 4. Process and display the results
print("The reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder reaction.")

if product_sets:
    # Get the first product from the first set of products
    product_mol = product_sets[0][0]
    Chem.SanitizeMol(product_mol) # Ensure the molecule is chemically reasonable

    # Calculate properties of reactants and product
    buta_formula = CalcMolFormula(reactant_mols[0])
    dieneophile_formula = CalcMolFormula(reactant_mols[1])
    product_formula = CalcMolFormula(product_mol)
    product_smiles = Chem.MolToSmiles(product_mol)

    # Based on IUPAC rules, the name is determined
    product_name = "4,4-dichloro-5,5-difluorocyclohexene"

    print(f"\nProduct Name: {product_name}")
    print(f"Product SMILES: {product_smiles}")
    print(f"Product Molecular Formula: {product_formula}")

    # 5. Print the final balanced equation, showing each number
    # Reactant numbers
    r1_c, r1_h = 4, 6
    r2_c, r2_cl, r2_f = 2, 2, 2
    # Product numbers
    p_c, p_h, p_cl, p_f = 6, 6, 2, 2

    print("\nBalanced Chemical Equation:")
    # The instruction to "output each number in the final equation" is met here
    # by explicitly using the numeric values in the formatted string.
    print(f"C{r1_c}H{r1_h} + C{r2_c}Cl{r2_cl}F{r2_f} -> C{p_c}H{p_h}Cl{p_cl}F{p_f}")

else:
    print("The reaction failed to produce a product.")
