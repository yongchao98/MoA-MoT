# This script explains the chemical reaction and identifies the products A and B.
# If you do not have the rdkit library installed, you can install it by running:
# pip install rdkit

try:
    from rdkit import Chem
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit'")
    print("Continuing without molecule visualization.\n")
    Chem = None

# --- Reaction Explanation ---
# The reaction is an iron-catalyzed radical addition of a peroxide to an alkene.

# 1. Initiation:
# The peroxide, tert-butyl peroxybenzoate, cleaves at its weak O-O bond,
# facilitated by heat and the Fe(III) catalyst. This forms a tert-butoxy
# radical (tBuO•) and a benzoyloxy radical (PhCOO•).

# 2. Propagation & Addition:
# These radicals add across the styrene double bond. The addition occurs at the
# terminal CH2 carbon to form the more stable benzylic radical. This leads
# to two possible radical intermediates.

# 3. Trapping & Product Formation:
# The benzylic radical intermediates are then 'trapped' by a functional group
# from the peroxide, leading to the final products. Since two pathways are
# possible, two major products (A and B), which are regioisomers, are formed.

# --- Product Identification ---

# Product A is formed when the tert-butoxy radical adds first, followed by trapping
# with a benzoyloxy group.
product_A_name = "2-(tert-butoxy)-1-phenylethyl benzoate"
product_A_smiles = "CC(C)(C)OCC(c1ccccc1)OC(=O)c2ccccc2"

# Product B is formed when the benzoyloxy radical adds first, followed by trapping
# with a tert-butoxy group.
product_B_name = "1-(tert-butoxy)-1-phenylethyl benzoate"
product_B_smiles = "CC(C)(C)OC(c1ccccc1)COC(=O)c2ccccc2"

# --- Output the results ---
print("The two major products, A and B, are the following regioisomers:\n")

print("Product A:")
print(f"  Name: {product_A_name}")
print(f"  SMILES notation: {product_A_smiles}")
if Chem:
    mol_A = Chem.MolFromSmiles(product_A_smiles)
    # The following line will display the structure in environments like Jupyter notebooks
    # In a terminal, it will just show that a molecule object was created.
    # display(mol_A)

print("\nProduct B:")
print(f"  Name: {product_B_name}")
print(f"  SMILES notation: {product_B_smiles}")
if Chem:
    mol_B = Chem.MolFromSmiles(product_B_smiles)
    # display(mol_B)
