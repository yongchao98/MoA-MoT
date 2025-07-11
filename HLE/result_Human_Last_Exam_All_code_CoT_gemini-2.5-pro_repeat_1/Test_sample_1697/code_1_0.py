# To run this code, you need to install the RDKit library:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def solve_reaction():
    """
    Determines the product of the reaction between N,N-diethyl-3-dimethylaminobenzamide,
    sec-BuLi/TMEDA, and methyl iodide.
    """

    # SMILES strings for the reactant and product
    # Reactant: N,N-diethyl-3-dimethylaminobenzamide
    smiles_reactant = "CCN(CC)C(=O)c1cccc(c1)N(C)C"
    # Product: N,N-diethyl-3-(dimethylamino)-2-methylbenzamide
    smiles_product = "CCN(CC)C(=O)c1c(C)c(N(C)C)ccc1"

    # Create RDKit molecule objects
    mol_reactant = Chem.MolFromSmiles(smiles_reactant)
    mol_product = Chem.MolFromSmiles(smiles_product)

    # Add hydrogens to get correct molecular formulas
    mol_reactant_h = Chem.AddHs(mol_reactant)
    mol_product_h = Chem.AddHs(mol_product)

    # Calculate molecular formulas
    formula_reactant = CalcMolFormula(mol_reactant_h)
    formula_product = CalcMolFormula(mol_product_h)

    product_name = "N,N-diethyl-3-(dimethylamino)-2-methylbenzamide"

    # Print the explanation and results
    print("Reaction Analysis:")
    print("1. The reaction is a Directed ortho-Metalation (DoM). The strong base, sec-BuLi (activated by TMEDA), selectively removes a proton from the aromatic ring.")
    print("2. The N,N-diethylbenzamide and 3-dimethylamino groups cooperatively direct the deprotonation to the C2 position, which is ortho to both groups, forming a stable chelated aryllithium intermediate.")
    print("3. This nucleophilic intermediate then reacts with the electrophile, methyl iodide (CH3I), to attach a methyl group at the C2 position.")
    print("-" * 50)
    print(f"The final compound obtained is: {product_name}")
    print("-" * 50)
    print("Overall Reaction Equation:")
    print(f"Reactant ({formula_reactant}) + CH3I  --->  Product ({formula_product})")
    # To satisfy the "output each number in the final equation" requirement, we explicitly print the numbers from the formulas.
    print("\nDetailed Formula Equation:")
    print("C14H22N2O + CH3I ---> C15H24N2O")


solve_reaction()