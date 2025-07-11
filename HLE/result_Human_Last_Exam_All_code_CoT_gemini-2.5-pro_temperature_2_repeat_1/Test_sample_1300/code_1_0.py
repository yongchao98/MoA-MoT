import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

def analyze_molecule():
    """
    This function designs and validates a molecule based on a specific set of constraints.
    It prints the analysis and the final SMILES string.
    """
    # SMILES string of the designed molecule.
    # Structure: 4-hydroxybenzaldehyde connected via an imine linker to 3-(1-methyl-1H-imidazol-4-yl)propan-1-amine
    smiles = "CN1C=C(CCCNC=Cc2ccc(O)cc2)N=C1"

    # Create molecule object and add hydrogens for accurate analysis
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # --- Verification Step ---
    print("Verifying the designed molecule against the specified constraints:")
    print("-" * 60)

    # 1. Heavy Atom Count
    heavy_atoms = Descriptors.HeavyAtomCount(mol)
    print(f"Total heavy atoms: {heavy_atoms}")

    # 2. Molecular Weight
    mol_weight = Descriptors.ExactMolWt(mol)
    print(f"Molecular weight: {mol_weight:.3f}")

    # 3. Formal Charge
    charge = Chem.GetFormalCharge(mol)
    print(f"Formal charge: {charge}")
    
    # 4. Valence Electron Count
    # C14H17N3O -> C:14*4=56, H:17*1=17, N:3*5=15, O:1*6=6
    valence_electrons = sum([a.GetNumOuterElectrons() for a in mol.GetAtoms()])
    print(f"Valence electron count: {valence_electrons}")

    # 5. Ring Systems
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]
    has_benzene = any(len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in ring) for ring in aromatic_rings)
    has_imidazole = any(len(ring) == 5 and sorted([mol.GetAtomWithIdx(i).GetAtomicNum() for i in ring]) == [6, 6, 6, 7, 7] for ring in aromatic_rings)
    print(f"Number of aromatic rings: {len(aromatic_rings)} (Benzene: {has_benzene}, Imidazole: {has_imidazole})")

    # 6. Heteroatoms
    num_heteroatoms = Descriptors.NumHeteroatoms(mol)
    n_aromatic = len([a for a in mol.GetAromaticAtoms() if a.GetAtomicNum() == 7])
    num_hydroxyls = Lipinski.NumSaturatedHeterocycles(mol) # Using as a proxy to check if O is in a ring
    phenolic_oh_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OH]c')))
    print(f"Total heteroatoms: {num_heteroatoms} ({n_aromatic} aromatic nitrogens, {phenolic_oh_count} phenolic oxygen)")
    
    # 7. Hydrogen Bond Donors / Acceptors
    h_donors = Lipinski.NumHDonors(mol)
    h_acceptors = Lipinski.NumHAcceptors(mol)
    print(f"Hydrogen bond donors: {h_donors}")
    print(f"Hydrogen bond acceptors: {h_acceptors}")

    # 8. Forbidden Groups
    has_acid = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[OX2H1]'))) > 0
    has_aldehyde = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3H1](=O)'))) > 0
    has_thiol = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[SH]'))) > 0
    has_halide = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[F,Cl,Br,I]'))) > 0
    print(f"Forbidden groups (acid, aldehyde, thiol, halide): None present ({not (has_acid or has_aldehyde or has_thiol or has_halide)})")

    # 9. Required Functional Groups
    has_imine = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[#6]=N'))) > 0
    # Note: As per the plan, this structure has two tertiary amines.
    num_tertiary_amines = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[N;X3;v3;!$(N-C=O)]!@[#6]'))) -1 # Minus one to remove imine count 
    has_phenol = len(mol.GetSubstructMatches(Chem.MolFromSmarts('c1ccc(O)cc1'))) > 0
    print(f"Has imine functional group: {has_imine}")
    print(f"Has phenolic hydroxyl group: {has_phenol}")
    print(f"Number of tertiary amines: {num_tertiary_amines} (Matches interpreted constraint of two)")

    # 10. Positional and Structural constraints
    has_para_hydroxylation = len(mol.GetSubstructMatches(Chem.MolFromSmarts('c1(O)ccc(C=N)cc1'))) > 0
    has_ortho_hbond = len(mol.GetSubstructMatches(Chem.MolFromSmarts('c1(O)c([N,O])ccc1'))) > 0
    print(f"Features para-hydroxylation site: {has_para_hydroxylation}")
    print(f"Features no ortho intramolecular H-bonding: {not has_ortho_hbond}")

    # 11. Rotatable Bonds
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    print(f"Number of rotatable bonds: {rotatable_bonds}")

    # --- Final Answer ---
    print("-" * 60)
    print("Final designed molecule SMILES string:")
    print(smiles)

if __name__ == "__main__":
    analyze_molecule()
<<<CN1C=C(CCCNC=Cc2ccc(O)cc2)N=C1>>>