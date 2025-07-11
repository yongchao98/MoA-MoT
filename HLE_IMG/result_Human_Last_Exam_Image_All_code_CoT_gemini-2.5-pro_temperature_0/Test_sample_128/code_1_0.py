# RDKit is a cheminformatics library. If you don't have it, install it using:
# pip install rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit'")
    print("Proceeding with pre-calculated values.")
    # Pre-calculated values for fallback
    product_name = "2-(pyridin-2-yl)isoindole-1-carbonitrile"
    molecular_formula = "C14H9N3"
    molecular_weight = 219.080
else:
    # The product is 2-(pyridin-2-yl)isoindole-1-carbonitrile.
    # We can represent this molecule using its SMILES string.
    product_smiles = "N#Cc1c2ccccc2n1c3ncccc3"
    mol_A = Chem.MolFromSmiles(product_smiles)
    
    product_name = "2-(pyridin-2-yl)isoindole-1-carbonitrile"
    molecular_formula = CalcMolFormula(mol_A)
    molecular_weight = Descriptors.ExactMolWt(mol_A)

# The reaction is a 1:1:1 condensation.
print("The chemical equation for the formation of Compound A is:")
print("1 (2-aminopyridine) + 1 (o-phthalaldehyde) + 1 (TMSCN) -> 1 (Compound A) + Byproducts")
print("\nReaction conditions:")
print(f"Temperature: 28 C")
print(f"Time: 4 h")
print("-" * 40)
print("The identity of the final product, Compound A, is:")
print(f"Name: {product_name}")
print(f"Molecular Formula: {molecular_formula}")
print(f"Exact Molecular Weight: {molecular_weight:.4f}")
print("-" * 40)