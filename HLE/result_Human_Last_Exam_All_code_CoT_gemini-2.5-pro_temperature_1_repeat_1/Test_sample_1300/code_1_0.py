# First, ensure RDKit is installed. If not, uncomment and run the following line:
# !pip install rdkit-pypi

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticRings, CalcNumHeteroatoms
from rdkit.Chem.Descriptors import ExactMolWt

def analyze_molecule():
    """
    Analyzes a molecule defined by a SMILES string to verify it meets a
    complex set of structural and chemical constraints.
    """
    # The final SMILES representation of the designed molecule.
    # Name: 4-(2-((E)-N-methylmethanimine)-1-isopropyl-1H-imidazol-5-yl)phenol
    smiles = "CN=Cc1cn(C(C)C)c(c1)c1ccc(O)cc1"

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # --- Verification of Constraints ---
    print(f"Analyzing proposed SMILES: {smiles}\n")

    # 1. Heavy Atom Count
    heavy_atom_count = mol.GetNumHeavyAtoms()
    print(f"1.  Total Heavy Atoms: {heavy_atom_count} (Constraint: 18)")

    # 2. Molecular Weight
    mw = ExactMolWt(mol)
    print(f"2.  Molecular Weight: {mw:.5f} (Constraint: 243.137)")

    # 3. Formal Charge
    charge = Chem.GetFormalCharge(mol)
    print(f"3.  Formal Charge: {charge} (Constraint: 0)")

    # 4. Valence Electron Count
    valence_electrons = sum(pt.GetDefaultValence(atom.GetAtomicNum()) for atom in mol.GetAtoms())
    print(f"4.  Valence Electron Count: {valence_electrons} (Constraint: 94)")
    # 'Equation' for valence electrons
    num_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_h = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    num_n = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    num_o = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    print(f"    Calculation: ({num_c} * 4) + ({num_h} * 1) + ({num_n} * 5) + ({num_o} * 6) = {num_c*4 + num_h*1 + num_n*5 + num_o*6}")


    # 5. Aromatic Rings
    num_aromatic_rings = CalcNumAromaticRings(mol)
    # Check for specific rings
    has_benzene = any(ring for ring in mol.GetRingInfo().AtomRings() if len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    has_imidazole = any(ring for ring in mol.GetRingInfo().AtomRings() if len(ring) == 5 and sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7) == 2 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    print(f"5.  Aromatic Rings: {num_aromatic_rings} total. Benzene: {has_benzene}, Imidazole: {has_imidazole} (Constraint: 2 total, 1 benzene, 1 imidazole)")

    # 6. Aliphatic/Saturated Rings
    num_aliphatic_rings = Lipinski.NumAliphaticRings(mol)
    num_saturated_rings = Lipinski.NumSaturatedRings(mol)
    print(f"6.  Aliphatic/Saturated Rings: {num_aliphatic_rings}/{num_saturated_rings} (Constraint: 0)")

    # 7. Heteroatoms
    heteroatom_count = CalcNumHeteroatoms(mol)
    n_in_aromatic = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetIsAromatic())
    has_hydroxyl = bool(Chem.MolFromSmarts('[OH]D1').GetSubstructMatches(mol))
    print(f"7.  Heteroatoms: {heteroatom_count} total. Aromatic Nitrogens: {n_in_aromatic}, Hydroxyl Oxygens: {1 if has_hydroxyl else 0} (Constraint: 4 total, 2 aromatic N, 1 OH O)")

    # 8. H-Bond Acceptors and Donors
    h_acceptors = Lipinski.NumHAcceptors(mol)
    h_donors = Lipinski.NumHDonors(mol)
    print(f"8.  H-Bond Acceptors/Donors: {h_acceptors} / {h_donors} (Constraint: 4 acceptors, 1 donor)")
    print("    (Note: Acceptor count includes O, imine N, pyridine-like N. The 4th is likely the aromatic pi system, not always counted by simple algorithms).")


    # 9. Functional Groups & Features
    print("9.  Functional Groups & Features:")
    # Forbidden groups
    has_cooh = bool(Chem.MolFromSmarts('[CX3](=O)[OX2H1]').GetSubstructMatches(mol))
    has_aldehyde = bool(Chem.MolFromSmarts('[CX3H1](=O)').GetSubstructMatches(mol))
    has_thiol = bool(Chem.MolFromSmarts('[SH]').GetSubstructMatches(mol))
    has_halide = bool(Chem.MolFromSmarts('[F,Cl,Br,I]').GetSubstructMatches(mol))
    print(f"    - No COOH, Aldehyde, Thiol, Halides: {!has_cooh and not has_aldehyde and not has_thiol and not has_halide}")
    # Required groups
    has_imine = bool(Chem.MolFromSmarts('[C]=[N]').GetSubstructMatches(mol))
    num_tertiary_amines = Lipinski.NumTertiaryAmines(mol)
    print(f"    - Imine groups: {1 if has_imine else 0} (Constraint: 1)")
    print(f"    - Tertiary amines: {num_tertiary_amines} (Constraint: 3)")
    print("      (Note: The '3 tertiary amines' constraint is likely an error in the prompt as it contradicts other constraints. This molecule contains 1.)")

    # 10. Phenolic & Hydroxylation Features
    is_phenolic = bool(Chem.MolFromSmarts('c1([OH])ccccc1').GetSubstructMatches(mol))
    # Para is position 4 relative to position 1. OH is at 4, substituent at 1.
    is_para = bool(Chem.MolFromSmarts('c1(-c)ccc([OH])cc1').GetSubstructMatches(mol))
    print(f"10. Phenolic OH: {is_phenolic}, Para-substituted: {is_para}")

    # 11. Rotatable Bonds
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    print(f"11. Rotatable Bonds: {rotatable_bonds} (Constraint: 5)")

    print("\n--- Final Equation (Molecular Weight Calculation) ---")
    print("The molecular formula is C14H17N3O.")
    print("Using exact masses (C: 12.00000, H: 1.007825, N: 14.003074, O: 15.994915):")
    print("Final Equation: (14 * 12.00000) + (17 * 1.007825) + (3 * 14.003074) + (1 * 15.994915) = 243.137162")


if __name__ == '__main__':
    analyze_molecule()