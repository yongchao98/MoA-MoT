# First, ensure you have RDKit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def solve_chemistry_problem():
    """
    This function solves the multi-step synthesis problem to identify compound C.
    """
    # Step 1: Define the structure of Compound B.
    # Based on the reaction sequence, compound B is (diethylamino)tris(2,4,6-trimethoxyphenyl)methane.
    # It has a central carbon bonded to a diethylamino group and three 2,4,6-trimethoxyphenyl groups.
    # Each trimethoxyphenyl group has 3 methoxy groups, so there are 3 * 3 = 9 methoxy groups in total.
    smiles_B = 'CCN(CC)C(c1c(OC)cc(OC)cc1OC)(c2c(OC)cc(OC)cc2OC)c3c(OC)cc(OC)cc3OC)'
    mol_B = Chem.MolFromSmiles(smiles_B)

    # Step 2: Define the demethylation reaction.
    # LiI cleaves aryl methyl ethers (Ar-O-CH3) to form phenols (Ar-OH).
    # The reaction SMARTS represents this transformation.
    rxn_smarts = '[c:1]-[O:2]-[CH3:3]>>[c:1]-[O:2][H]'
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)

    # Step 3: Apply the reaction to Compound B to get Compound C.
    # We loop the reaction to ensure all 9 methoxy groups are demethylated.
    product_mol_C = mol_B
    while True:
        products = rxn.RunReactants((product_mol_C,))
        if not products:
            break
        # Take the first product from the reaction output
        product_mol_C = products[0][0]

    # Sanitize the final molecule to ensure valences are correct.
    Chem.SanitizeMol(product_mol_C)

    # Step 4: Output the results for Compound C.
    smiles_C = Chem.MolToSmiles(product_mol_C)
    
    # Calculate molecular formula
    mol_C_with_H = Chem.AddHs(product_mol_C)
    formula_C = Descriptors.rdMolDescriptors.CalcMolFormula(mol_C_with_H)

    print("The final product, Compound C, is (diethylamino)tris(2,4,6-trihydroxyphenyl)methane.")
    print(f"SMILES string of Compound C: {smiles_C}")
    print(f"Molecular Formula of Compound C: {formula_C}")

    # Outputting each number in the final formula equation
    atom_counts = Descriptors.rdMolDescriptors.CalcComposition(mol_C_with_H)
    print("\nElemental Composition of Compound C:")
    for atom_symbol, count in atom_counts:
        # The prompt requires outputting each number in the final equation.
        # This will print the count for each element.
        print(f"Number of {Chem.Atom(int(atom_symbol)).GetSymbol()} atoms: {int(count)}")

if __name__ == '__main__':
    solve_chemistry_problem()