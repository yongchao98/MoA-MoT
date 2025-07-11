import sys

# This script uses the RDKit library.
# If you do not have it installed, you can install it via pip:
# pip install rdkit-pypi

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("Error: RDKit library not found.", file=sys.stderr)
    print("Please install it using the command: pip install rdkit-pypi", file=sys.stderr)
    sys.exit(1)

# Step-by-step reaction simulation
# Starting Material: Terpinolene
# SMILES representation of 1-methyl-4-(propan-2-ylidene)cyclohex-1-ene
terpinolene_smiles = 'CC1=CCC(=C(C)C)CC1'
terpinolene = Chem.MolFromSmiles(terpinolene_smiles)

# Step 1: Epoxidation of the endocyclic double bond with m-CPBA
# We use a reaction SMARTS to represent this transformation
# [C;R:1]=[C;R:2] means a double bond between two atoms in a ring
# >>[C:1]1O[C:2]1 means those two atoms form an epoxide ring
rxn_step1 = AllChem.ReactionFromSmarts('[C;R:1]=[C;R:2]>>[C:1]1O[C:2]1')
products_step1 = rxn_step1.RunReactants((terpinolene,))
compound1 = products_step1[0][0]
Chem.SanitizeMol(compound1)
compound1_smiles = Chem.MolToSmiles(compound1)
print(f"Step 1 Product (Compound 1): 1,2-epoxy-p-menth-4(8)-ene\nSMILES: {compound1_smiles}\n")


# Step 2: Conversion of epoxide to thiirane with N,N-dimethyl thioformamide
# Reaction SMARTS to replace the oxygen of the epoxide with sulfur
rxn_step2 = AllChem.ReactionFromSmarts('[C:1]1O[C:2]1>>[C:1]1S[C:2]1')
products_step2 = rxn_step2.RunReactants((compound1,))
compound2 = products_step2[0][0]
Chem.SanitizeMol(compound2)
compound2_smiles = Chem.MolToSmiles(compound2)
print(f"Step 2 Product (Compound 2): 1,2-epithio-p-menth-4(8)-ene\nSMILES: {compound2_smiles}\n")

# Step 3: Reduction of thiirane with LiAlH4
# This opens the thiirane ring. Hydride attacks the less hindered carbon.
# The result is a thiol at the more substituted carbon.
# The SMARTS pattern targets a thiirane ring with a degree-4 and degree-3 carbon
rxn_step3 = AllChem.ReactionFromSmarts('[C;D4:1]1S[C;D3:2]1>>[C:1](S)[C:2]')
products_step3 = rxn_step3.RunReactants((compound2,))
compound3 = products_step3[0][0]
Chem.SanitizeMol(compound3)
compound3_smiles = Chem.MolToSmiles(compound3, canonical=True)
compound3_name = "1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol"

print(f"Step 3 Product (Compound 3): {compound3_name}\nSMILES: {compound3_smiles}")
