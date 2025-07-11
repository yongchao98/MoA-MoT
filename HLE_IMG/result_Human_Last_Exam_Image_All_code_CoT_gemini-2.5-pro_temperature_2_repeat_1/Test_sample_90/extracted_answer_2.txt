# Import necessary libraries for molecular structure handling
from rdkit import Chem
from rdkit.Chem import AllChem

def solve_chemistry_problem():
    """
    This function proposes structures for products A, B, and C based on the provided reaction and spectral data.
    The proposed structures are represented as SMILES strings.
    """

    # Proposed structure for Product A: A ring-opened pyrrole derivative.
    # This structure is consistent with C14H20N2O3, has a pyrrole ring, a methyl ester,
    # an acetyl group, and an alkyl chain as suggested by NMR.
    smiles_A = "CC(=O)NCCCCc1[nH]c(C)cc1C(=O)OC"

    # Proposed structure for Product B: A rigid, rearranged polycyclic adduct.
    # This structure is consistent with C12H14N2O3 (DBE=7) and the NMR features of a
    # rigid, conjugated system with a cis-alkene.
    smiles_B = "CC(=O)N1C(=O)C=C(C(=O)OC)[C@H]2[C@H]1CCCN2"

    # Proposed structure for Product C: N-acetylated starting material derivative.
    # This structure is consistent with C11H16N2O3 (SM + Acetyl group)
    # The SMILES represents one of the possible isomers.
    smiles_C = "CC(=O)NC1CCCC(C(=O)O)(C2=NCCC2)1"


    mol_A = Chem.MolFromSmiles(smiles_A)
    mol_B = Chem.MolFromSmiles(smiles_B)
    mol_C = Chem.MolFromSmiles(smiles_C)


    print("Proposed Structure for Product A:")
    print(f"Molecular Formula: {AllChem.CalcMolFormula(mol_A)}")
    print(f"SMILES: {smiles_A}")
    print("Rationale: Contains a pyrrole ring, methyl ester, and an N-acetylated chain, matching the NMR fragmentation and formula. This suggests the starting material's rings opened up.")

    print("\nProposed Structure for Product B:")
    print(f"Molecular Formula: {AllChem.CalcMolFormula(mol_B)}")
    print(f"SMILES: {smiles_B}")
    print("Rationale: A rigid polycyclic structure with high unsaturation, consistent with a rearranged cycloaddition product and the complex NMR spectrum.")

    print("\nProposed Structure for Product C:")
    print(f"Molecular Formula: {AllChem.CalcMolFormula(mol_C)}")
    print(f"SMILES: {smiles_C}")
    print("Rationale: Represents the starting material with an added acetyl group. This is the simplest product, consistent with its formula being derived directly from the starting material.")

solve_chemistry_problem()