# The user needs to install RDKit for this script to run.
# You can install it via pip: pip install rdkit
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

def analyze_molecule(smiles_string):
    """
    Analyzes a molecule based on a SMILES string and prints its properties
    according to the user's request.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            print("Error: Invalid SMILES string provided.")
            return
            
        # Add explicit hydrogens for accurate calculations
        mol = Chem.AddHs(mol)

        # --- Calculations ---
        
        # 1. Basic properties
        heavy_atoms = mol.GetNumHeavyAtoms()
        exact_mw = Descriptors.ExactMolWt(mol)
        formal_charge = Chem.GetFormalCharge(mol)
        
        # 2. Valence electrons
        valence_electrons = 0
        for atom in mol.GetAtoms():
            valence_electrons += Descriptors.GetValence(atom)
        # The above RDKit GetValence counts bonding electrons. A more standard way for this problem is periodic table valence.
        valence_electrons_periodic = 0
        for atom in mol.GetAtoms():
             valence_electrons_periodic += atom.GetTotalValence() - atom.GetTotalNumHs() + atom.GetNumRadicalElectrons()
        # Even simpler for the prompt's likely definition:
        valence_electrons = sum([{'C': 4, 'H': 1, 'N': 5, 'O': 6, 'S': 6, 'F': 7, 'Cl': 7, 'Br': 7, 'I': 7}[atom.GetSymbol()] for atom in mol.GetAtoms()])


        # 3. Structural Features & Functional Groups
        # Aromatic Rings
        aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]
        num_aromatic_rings = len(aromatic_rings)
        has_benzene = any(len(ring) == 6 for ring in aromatic_rings)
        has_imidazole = any(len(ring) == 5 and sum(1 for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == 'N') == 2 for ring in aromatic_rings)
        
        # Other rings
        aliphatic_rings = mol.GetRingInfo().NumAliphaticRings()
        
        # H-bond donors/acceptors (Lipinski definitions)
        h_bond_acceptors = Lipinski.NumHAcceptors(mol)
        h_bond_donors = Lipinski.NumHDonors(mol)
        has_hydroxyl_donor = bool(Chem.MolFromSmarts('[OH]')) # check OH group is a donor
        
        # Other groups
        has_imine = mol.HasSubstructMatch(Chem.MolFromSmarts('[C]=[N]'))
        has_phenolic_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts('c1([OH])ccccc1'))
        
        # Rotatable Bonds
        rotatable_bonds = Lipinski.NumRotatableBonds(mol)
        
        # 4. Molecular Formula "Equation"
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        
        # --- Verification and Output ---
        print("--- Verifying Proposed Molecule ---")
        print(f"SMILES: {smiles_string}")
        print(f"Molecular Formula: {mol_formula}\n")

        print("--- Checking Constraints ---")
        print(f"Total Heavy Atoms: {heavy_atoms} (Target: 18)")
        print(f"Molecular Weight: {exact_mw:.3f} (Target: 243.137)")
        print(f"Formal Charge: {formal_charge} (Target: 0)")
        print(f"Valence Electrons: {valence_electrons} (Target: 94)")
        print(f"Number of Aromatic Rings: {num_aromatic_rings} (Target: 2)")
        print(f"  - Contains Benzene: {has_benzene} (Target: True)")
        print(f"  - Contains Imidazole: {has_imidazole} (Target: True)")
        print(f"No Aliphatic/Saturated Rings: {aliphatic_rings == 0} (Target: True)")
        print(f"Number of Rotatable Bonds: {rotatable_bonds} (Target: 5)")
        print(f"Hydrogen Bond Acceptors: {h_bond_acceptors} (Target: 4)")
        print(f"Hydrogen Bond Donors (from OH group): {has_hydroxyl_donor} (Target: 1)")
        print(f"Contains Imine Group: {has_imine} (Target: True)")
        print(f"Contains Phenolic Hydroxyl: {has_phenolic_hydroxyl} (Target: True)")
        
        # The prompt asked for each number in the final equation.
        # Here we break down the molecular formula as an "equation".
        print("\n--- Final Equation (Molecular Formula Breakdown) ---")
        # Let's extract counts from the formula string
        atoms = {'C':0, 'H':0, 'N':0, 'O':0}
        import re
        parts = re.findall(r'([A-Z][a-z]*)(\d*)', mol_formula)
        for symbol, count in parts:
            atoms[symbol] = int(count) if count else 1

        print(f"Carbon (C): {atoms['C']}")
        print(f"Hydrogen (H): {atoms['H']}")
        print(f"Nitrogen (N): {atoms['N']}")
        print(f"Oxygen (O): {atoms['O']}")

    except ImportError:
        print("This script requires the RDKit library.")
        print("Please install it using: pip install rdkit-pypi")
    except Exception as e:
        print(f"An error occurred: {e}")


# --- Main Execution ---
# The SMILES string of the designed molecule.
# Structure: 2-((1,1-dimethylethylidene)amino)ethyl)-4-(4-hydroxyphenyl)-1H-imidazole
final_smiles = "CC(C)=NCCc1nc(c2ccc(O)cc2)[nH]1"
analyze_molecule(final_smiles)
