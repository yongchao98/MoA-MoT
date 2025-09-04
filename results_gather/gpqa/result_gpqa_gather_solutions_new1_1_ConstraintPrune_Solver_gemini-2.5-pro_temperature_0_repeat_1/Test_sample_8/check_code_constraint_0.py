import sys
import subprocess

# --- Installation Check ---
# This part ensures the necessary libraries are installed.
try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("RDKit is not installed. Attempting to install...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
        print("RDKit installed successfully. Please re-run the script.")
    except Exception as e:
        print(f"Failed to install RDKit: {e}")
    sys.exit()

try:
    import pyscf
    from pyscf import gto, symm
except ImportError:
    print("PySCF is not installed. Attempting to install...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyscf"])
        print("PySCF installed successfully. Please re-run the script.")
    except Exception as e:
        print(f"Failed to install PySCF: {e}")
    sys.exit()

import numpy as np

def get_point_group(smiles, name, num_confs=1, tol=1e-2):
    """
    Generates a 3D structure from SMILES and determines its point group using PySCF.
    For flexible molecules, it checks multiple conformers to find the highest symmetry.
    """
    try:
        mol_rdkit = Chem.MolFromSmiles(smiles)
        if not mol_rdkit:
            return f"Error: RDKit could not parse SMILES for {name}"
        mol_rdkit = Chem.AddHs(mol_rdkit)

        # Generate conformer(s)
        if num_confs > 1:
            cids = AllChem.EmbedMultipleConfs(mol_rdkit, numConfs=num_confs, params=AllChem.ETKDGv3())
            if cids:
                AllChem.MMFFOptimizeMoleculeConfs(mol_rdkit)
        else:
            cid = AllChem.EmbedMolecule(mol_rdkit, AllChem.ETKDGv3())
            if cid != -1:
                AllChem.MMFFOptimizeMolecule(mol_rdkit)

        if mol_rdkit.GetNumConformers() == 0:
            return f"Error: RDKit could not generate a conformer for {name}"

        highest_symmetry_group = 'C1'
        highest_symmetry_order = 1

        for conf_id in range(mol_rdkit.GetNumConformers()):
            conf = mol_rdkit.GetConformer(conf_id)
            atoms_pyscf = []
            for i in range(mol_rdkit.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                atom_symbol = mol_rdkit.GetAtomWithIdx(i).GetSymbol()
                atoms_pyscf.append([atom_symbol, (pos.x, pos.y, pos.z)])

            mol_pyscf = gto.Mole()
            mol_pyscf.atom = atoms_pyscf
            mol_pyscf.build(0, 0)
            
            pg = symm.topgroup(mol_pyscf, tol=tol)
            
            current_order = symm.order(pg) if pg in symm.GROUP_ORDER else 1
            if current_order > highest_symmetry_order:
                highest_symmetry_group = pg
                highest_symmetry_order = current_order
                
        return highest_symmetry_group

    except Exception as e:
        return f"Analysis Error for {name}: {e}"

def check_answer():
    """
    Checks the correctness of the LLM's answer by computationally determining
    the point group of each molecule.
    """
    llm_answer = 'C'

    molecules_to_check = {
        "A": ("triisopropyl borate", "CC(C)OB(OC(C)C)OC(C)C", 20, 5e-2),
        "B": ("benzo[...]trifuran[...]hexaone", "C1=C2C(=O)OC(=O)C2=C3C4=C(C(=C1)C(=O)OC4=O)C(=O)OC3=O", 1, 1e-2),
        "C": ("triphenyleno[...]trifuran[...]hexaone", "C1=C2C(=O)OC(=O)C2=C3C4=C(C=C5C(=O)OC(=O)C5=C4)C6=C3C=C1C(=O)OC6=O", 1, 1e-2),
        "D": ("quinuclidine", "C1CN2CCC1CC2", 1, 1e-2)
    }

    print("--- Starting Computational Verification ---")
    results = {}
    has_errors = False
    for key, (name, smiles, confs, tol) in molecules_to_check.items():
        print(f"Analyzing molecule {key}: {name}...")
        pg = get_point_group(smiles, name, num_confs=confs, tol=tol)
        if "Error" in str(pg):
            has_errors = True
        results[key] = pg
        print(f"Result for {key}: Point group found = {pg}")
    
    print("\n--- Analysis Summary ---")
    for key, pg in results.items():
        print(f"Molecule {key} ({molecules_to_check[key][0]}): {pg}")
    
    if has_errors:
        print("\n--- Verdict ---")
        print("Could not complete verification due to errors in structure generation or analysis.")
        return

    # --- Verdict Logic ---
    is_c_c3h = results.get('C') == 'C3h'
    is_b_d3h = results.get('B') == 'D3h'
    is_d_c3v = results.get('D') == 'C3v'
    is_a_not_c3h = results.get('A') != 'C3h'

    print("\n--- Verdict ---")
    if llm_answer != 'C':
        print(f"Incorrect. The provided answer is '{llm_answer}', but the analysis points to 'C' as the correct molecule with C3h symmetry.")
        return

    if is_c_c3h and is_b_d3h and is_d_c3v and is_a_not_c3h:
        print("Correct")
    else:
        reason = f"Incorrect. The provided answer is C, but the computational analysis yields conflicting results.\n"
        if not is_c_c3h:
            reason += f"- The target molecule C was found to have {results.get('C')} symmetry, not C3h as expected.\n"
        if not is_b_d3h:
            reason += f"- Molecule B (mellitic trianhydride), a known D3h molecule, was found to have {results.get('B')} symmetry.\n"
        if not is_d_c3v:
            reason += f"- Molecule D (quinuclidine), a known C3v molecule, was found to have {results.get('D')} symmetry.\n"
        if not is_a_not_c3h:
            reason += f"- Molecule A (triisopropyl borate) was unexpectedly found to have C3h symmetry in its lowest energy conformer.\n"
        
        # Check if another molecule was found to be C3h instead of C
        other_c3h = [k for k, v in results.items() if v == 'C3h' and k != 'C']
        if other_c3h:
            reason += f"The analysis suggests that molecule(s) {', '.join(other_c3h)} are the actual C3h species."
        
        print(reason.strip())

# Run the check
check_answer()