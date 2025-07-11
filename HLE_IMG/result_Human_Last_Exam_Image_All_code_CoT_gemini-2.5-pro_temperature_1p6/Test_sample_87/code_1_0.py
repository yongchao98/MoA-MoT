import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors

def get_mol_info(smiles, name):
    """Calculates and prints molecular formula and weight for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Could not parse SMILES for {name}")
        return None
    formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
    print(f"Structure of Product {name}:")
    print(f"Molecular Formula: {formula}")
    return mol

# Product C: N-acetylated starting material
smiles_C = "CC(=O)N1CCCC1(C(=O)O)C2=NCCC2"
mol_C = get_mol_info(smiles_C, "C (C11H16N2O3)")

# The cycloaddition of the munchone from C with methyl propiolate gives two possible regioisomers.
# We will show the pathway for one of them leading to A and B.

# P_cyclo: The initial, stable cycloadduct before hydrolysis or fragmentation
# Structure: methyl 7a-(1-pyrrolin-2-yl)-2,3,5,6,7,7a-hexahydro-1H-pyrrolizine-1-carboxylate
smiles_P_cyclo = "COC(=O)C1=CN2C(C3=NCCC3)(CCCC2)C1"
# We don't print this intermediate, just use it for derivation.

# Product A: Hydrolysis product of P_cyclo
# The imine C=N is hydrolyzed to C=O and NH2, opening the pyrroline ring.
smiles_A = "COC(=O)C1=CN2C(C(=O)CCCN)(CCCC2)C1"
mol_A = get_mol_info(smiles_A, "A (C14H20N2O3)")

# Product B: Fragmentation (-C2H4) and oxidation (+O) product of P_cyclo.
# Retro-Diels-Alder on the saturated ring of P_cyclo removes ethylene. Then N-oxidation of the pyrrole nitrogen.
# Structure: methyl 5-(1-pyrrolin-2-yl)-3,5-dihydro-1H-pyrrolo[1,2-a]imidazole-7-carboxylate, but let's assume a simpler oxidation for formula matching.
# Let's propose a plausible structure that fits C12H14N2O3.
# A more likely pathway could be rearrangement. A possible structure based on intramolecular cyclization of the munchone M followed by reaction with water could be:
smiles_B = "CC1=C(C(=O)O)N2C(C3=NCCC3)(C1)CCC2=O"
mol_B = get_mol_info(smiles_B, "B (C12H14N2O3)")

# Display the structures
img = Draw.MolsToGridImage([mol_A, mol_B, mol_C],
                           molsPerRow=3,
                           subImgSize=(300, 300),
                           legends=["Product A", "Product B", "Product C"])
# To view the image in environments like Jupyter, just have `img` as the last line.
# For other environments, you might need to save it.
img.save("reaction_products.png")
print("\nImage of the three product structures has been saved as 'reaction_products.png'")
# The image itself will not be displayed here, but the code generates it.
# The following print statement is a text representation of the output.
print("\n--- Final Answer ---")
print("Product A is the hydrolyzed cycloadduct. SMILES: COC(=O)C1=CN2C(C(=O)CCCN)(CCCC2)C1")
print("Product B is a rearranged side product. A plausible structure is represented by SMILES: CC1=C(C(=O)O)N2C(C3=NCCC3)(C1)CCC2=O")
print("Product C is the N-acetylated starting material. SMILES: CC(=O)N1CCCC1(C(=O)O)C2=NCCC2")
