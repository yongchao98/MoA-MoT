import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def install(package):
    """Installs a package using pip."""
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Ensure rdkit is installed
try:
    from rdkit import Chem
except ImportError:
    print("RDKit not found. Installing...")
    install("rdkit-pypi")
    from rdkit import Chem
    print("RDKit installed successfully.")

def check_answer():
    """
    Checks the correctness of the proposed answer by verifying its properties
    against the given spectral data.
    """
    # --- 1. Define Given Data and Candidate Structures ---
    
    # Spectral data from the question
    ms_m_plus = 156
    ms_m_plus_2_ratio = 0.32  # 32%
    ir_oh_broad = (2700, 3500)
    ir_co_sharp = 1720
    nmr_acid_proton = 11.0
    nmr_aromatic_signals = [
        {'ppm': 8.02, 'type': 'd', 'integration': 2},
        {'ppm': 7.72, 'type': 'd', 'integration': 2}
    ]

    # Candidate molecules
    candidates = {
        'A': {'name': 'Phenyl chloroformate', 'smiles': 'O=C(Cl)OC1=CC=CC=C1'},
        'B': {'name': '4-chlorobenzoic acid', 'smiles': 'O=C(O)c1ccc(Cl)cc1'},
        'C': {'name': '3-Chloro-2-hydroxybenzaldehyde', 'smiles': 'O=Cc1cccc(Cl)c1O'},
        'D': {'name': '2-chlorobenzoic acid', 'smiles': 'O=C(O)c1ccccc1Cl'}
    }
    
    # The proposed answer to check
    proposed_answer_key = 'B'
    proposed_molecule_info = candidates[proposed_answer_key]
    mol = Chem.MolFromSmiles(proposed_molecule_info['smiles'])
    
    # List to store reasons for incorrectness
    errors = []

    # --- 2. Mass Spectrometry Check ---
    # Check for the presence of one Chlorine atom
    formula = CalcMolFormula(mol)
    if formula.count('Cl') != 1:
        errors.append(f"MS Check Failed: Expected 1 Chlorine atom for the M+2 peak, but found {formula.count('Cl')} in {proposed_molecule_info['name']}.")
    
    # Check molecular weight for the most abundant isotopes (C=12, H=1, O=16, Cl=35)
    atomic_weights = {'C': 12, 'H': 1, 'O': 16, 'Cl': 35}
    mol_weight = 0
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in atomic_weights:
            mol_weight += atomic_weights[symbol]
    
    # The molecular weight is calculated based on integer masses of the most abundant isotopes.
    # We check if the calculated weight matches the M+ peak.
    if mol_weight != ms_m_plus:
        errors.append(f"MS Check Failed: The M+ peak is at m/z {ms_m_plus}, but the calculated molecular weight for {proposed_molecule_info['name']} (with 35Cl) is {mol_weight}.")

    # --- 3. IR Spectroscopy Check ---
    # Check for a carboxylic acid group, which causes the characteristic broad O-H and C=O peaks.
    carboxylic_acid_smarts = '[CX3](=O)[OX2H1]'
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(carboxylic_acid_smarts)):
        errors.append(f"IR Check Failed: The spectrum indicates a carboxylic acid (broad peak 3500-2700 cm-1, sharp peak ~1720 cm-1), but {proposed_molecule_info['name']} does not contain a -COOH group.")

    # --- 4. 1H NMR Spectroscopy Check ---
    # Check for the carboxylic acid proton
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(carboxylic_acid_smarts)):
        errors.append(f"NMR Check Failed: The spectrum shows a peak at {nmr_acid_proton} ppm (carboxylic acid proton), but the structure lacks this functional group.")

    # Check the aromatic region pattern
    aromatic_protons = []
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() and atom.GetAtomicNum() == 6: # Aromatic Carbon
            # Count attached hydrogens
            hydrogens = atom.GetTotalNumHs()
            if hydrogens > 0:
                aromatic_protons.append(atom)
    
    total_aromatic_protons = sum(a.GetTotalNumHs() for a in aromatic_protons)
    expected_aromatic_protons = sum(sig['integration'] for sig in nmr_aromatic_signals)

    if total_aromatic_protons != expected_aromatic_protons:
        errors.append(f"NMR Check Failed: The spectrum shows {expected_aromatic_protons} aromatic protons, but the structure has {total_aromatic_protons}.")
    else:
        # Check for the specific para-substitution pattern (2 sets of 2 equivalent protons)
        # We can do this by checking the symmetry of the aromatic protons
        # GetSymmetryClasses returns a tuple of atom indices for each symmetry class
        symmetry_classes = Chem.GetSymmSSSR(mol)
        
        # For 4-chlorobenzoic acid, we expect the 4 aromatic protons to fall into 2 symmetry classes of 2 protons each.
        # This leads to the observed "2 doublets, 2H each" pattern.
        # For 2-chlorobenzoic acid, all 4 aromatic protons would be unique (4 classes of 1 proton).
        
        # A simpler check for para-substitution:
        substituents = []
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 1:
            ring_atoms = list(ring_info.AtomRings()[0])
            for atom_idx in ring_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in ring_atoms: # This is a substituent
                        substituents.append(atom_idx)
            
            # In a 6-membered ring, para substituents are at positions i and i+3
            if len(substituents) == 2:
                pos1 = ring_atoms.index(substituents[0])
                pos2 = ring_atoms.index(substituents[1])
                if abs(pos1 - pos2) != 3:
                    errors.append(f"NMR Check Failed: The aromatic pattern (2 doublets, 2H each) indicates para-substitution, but the structure {proposed_molecule_info['name']} is not para-substituted.")
            else:
                 errors.append(f"NMR Check Failed: Expected a di-substituted ring, but found {len(substituents)} substituents.")
        else:
            errors.append("NMR Check Failed: Could not analyze ring structure.")


    # --- 5. Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Run the check and print the result
result = check_answer()
print(result)