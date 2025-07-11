import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
except ImportError:
    print("RDKit is not installed. Please install it using: pip install rdkit-pypi")
    sys.exit(1)

def analyze_molecule(smiles_string):
    """
    Analyzes a molecule based on a SMILES string and prints its properties
    against a predefined set of constraints.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print(f"Error: Could not parse SMILES string: {smiles_string}")
        return

    # Add hydrogens to the molecule to get correct properties
    mol = Chem.AddHs(mol)

    print(f"--- Analysis for SMILES: {smiles_string} ---\n")

    # --- Molecular Formula and Weight ---
    heavy_atoms = Descriptors.HeavyAtomCount(mol)
    mw = Descriptors.ExactMolWt(mol)
    charge = Chem.GetFormalCharge(mol)
    
    # Calculate valence electrons
    valence_electrons = 0
    for atom in mol.GetAtoms():
        valence_electrons += Descriptors.calcImplicitValence(atom)
        valence_electrons += atom.GetTotalNumHs()
        
    print("--- Core Properties ---")
    print(f"Total Heavy Atoms:       {heavy_atoms} (Target: 18)")
    print(f"Molecular Weight:        {mw:.5f} (Target: 243.137)")
    print(f"Formal Charge:           {charge} (Target: 0)")
    print(f"Valence Electron Count:  {valence_electrons} (Target: 94)")
    
    # --- Molecular Weight Calculation Breakdown ---
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mol_formula = Chem.MolToMolBlock(mol)
    counts = { 'C': 0, 'H': 0, 'N': 0, 'O': 0 }
    for line in mol_formula.split('\n'):
        if 'M  END' in line:
            break
        parts = line.split()
        if len(parts) > 4:
            atom_symbol = parts[3]
            if atom_symbol in counts:
                counts[atom_symbol] += 1
    
    print("\n--- Exact Mass Calculation ---")
    print(f"Molecular Formula: {formula}")
    print(f"({counts['C']} * 12.00000) + ({counts['H']} * 1.007825) + ({counts['N']} * 14.003074) + ({counts['O']} * 15.994915) = {mw:.5f}")


    # --- Structural Features ---
    ring_info = mol.GetRingInfo()
    num_aromatic_rings = sum(1 for ring in ring_info.AtomRings() if Chem.IsAromatic(mol, ring))
    
    heteroatoms = [a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() not in [1, 6]]
    aromatic_nitrogens = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[n]')))
    
    print("\n--- Structural Features ---")
    print(f"Aromatic Rings:          {num_aromatic_rings} (Target: 2)")
    print(f"Total Heteroatoms:       {len(heteroatoms)} (Target: 4)")
    print(f"Aromatic Nitrogens:      {aromatic_nitrogens} (Target: 2)")


    # --- Functional Groups and H-Bonding ---
    hbd = Lipinski.NumHBD(mol)
    hba = Lipinski.NumHBA(mol)
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    
    # SMARTS patterns for checking groups
    imine_pattern = Chem.MolFromSmarts('[CX3]=[NX2]')
    phenol_pattern = Chem.MolFromSmarts('c[OX2H]')
    tert_amine_pattern = Chem.MolFromSmarts('[NX3;H0;!$(*=[O,S,P])](C)(C)C') # classic tertiary amine
    
    has_imine = "Yes" if mol.HasSubstructMatch(imine_pattern) else "No"
    has_phenol = "Yes" if mol.HasSubstructMatch(phenol_pattern) else "No"
    num_tert_amines = len(mol.GetSubstructMatches(tert_amine_pattern))

    print("\n--- Functional Groups & Bonding ---")
    print(f"Hydrogen Bond Donors:    {hbd} (Target: 1 specified as OH)")
    print(f"Hydrogen Bond Acceptors: {hba} (Target: 4)")
    print(f"Rotatable Bonds:         {rotatable_bonds} (Target: 5)")
    print(f"Contains Imine Group:    {has_imine} (Target: Yes)")
    print(f"Contains Phenolic OH:    {has_phenol} (Target: Yes)")
    print(f"Tertiary Amines:         {num_tert_amines} (Target: 3 - NOTE: This constraint is likely erroneous)")

    # --- Forbidden Group Check ---
    forbidden = {
        "Carboxylic Acid": "[CX3](=O)[OX2H1]",
        "Aldehyde": "[CX3H1](=O)[#6]",
        "Thiol": "[SH]",
        "Halide": "[#9,#17,#35,#53]"
    }
    found_forbidden = "None"
    for name, smarts in forbidden.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            found_forbidden = name
            break
            
    print(f"Forbidden Groups Check:  {found_forbidden} (Target: None)")
    print("\n------------------------------------")


# This SMILES string represents the molecule designed to fit the constraints as closely as possible.
# Structure: A phenol ring is connected via a para-positioned -C=N-CH2- linker to an
# imidazole ring which is substituted with a propyl group.
final_smiles = "Oc1ccc(C=NCc2c(CCC)n[nH]c2)cc1"

# Run the analysis
analyze_molecule(final_smiles)

# The final answer as requested by the user format.
print("\nFinal proposed SMILES string:")
print(f"<<<{final_smiles}>>>")
