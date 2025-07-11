import sys
import subprocess

# This task requires rdkit, networkx, and numpy.
# We will check for their installation and install them if they are missing.
def install_and_import(package, import_name=None):
    if import_name is None:
        import_name = package
    try:
        __import__(import_name)
    except ImportError:
        print(f"Installing {package}... please wait.")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        except Exception as e:
            print(f"Error installing {package}. Please install it manually using 'pip install {package}'")
            print(f"Error details: {e}")
            sys.exit(1)
    finally:
        globals()[import_name] = __import__(import_name)

def solve_cheminformatics_problem():
    """
    Solves the multi-step cheminformatics problem as described.
    """
    # Ensure required packages are installed
    install_and_import('rdkit', 'rdkit')
    install_and_import('numpy', 'numpy')
    install_and_import('networkx', 'networkx')

    from rdkit import Chem
    from rdkit.Chem import GraphDescriptors
    from networkx.algorithms.polynomials import matching_polynomial

    # --- Step 1: Identify the reference BCKDH complex substrate ---
    print("--- Step 1: Identifying the Reference BCKDH Substrate ---")

    bckdh_substrates = {
        'KIV': 'CC(C)C(=O)C(=O)O',    # a-Ketoisovalerate
        'KIC': 'CC(C)CC(=O)C(=O)O',  # a-Ketoisocaproate
        'KMV': 'CCC(C)C(=O)C(=O)O'     # a-Keto-b-methylvalerate
    }

    bckdh_bertz = {}
    print("Calculating Bertz's complexity for BCKDH substrates...")
    for name, smi in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smi)
        bertz_ct = GraphDescriptors.BertzCT(mol)
        bckdh_bertz[name] = bertz_ct
        print(f"  - Bertz Complexity for {name}: {bertz_ct:.4f}")

    bertz_values = sorted(bckdh_bertz.values())
    median_bertz = numpy.median(bertz_values)
    print(f"\nSorted Bertz values: {[f'{v:.4f}' for v in bertz_values]}")
    print(f"The median Bertz complexity is: {median_bertz:.4f}")

    # Find the molecule(s) with the median value.
    # It might be one or more. We'll proceed by calculating Balaban J for it/them.
    ref_molecule_smi = ""
    for name, val in bckdh_bertz.items():
        if val == median_bertz:
            ref_molecule_name = name
            ref_molecule_smi = bckdh_substrates[name]
            # Since KIC and KMV have the same value, we just need one to proceed.
            break

    print(f"The BCKDH substrate with median complexity is '{ref_molecule_name}'.")

    # --- Step 2: Identify the target substance among nucleosides ---
    print("\n--- Step 2: Identifying the Target Nucleoside ---")
    ref_mol = Chem.MolFromSmiles(ref_molecule_smi)
    ref_balaban_j = GraphDescriptors.BalabanJ(ref_mol)
    print(f"Calculating Balaban J index for reference substrate '{ref_molecule_name}': {ref_balaban_j:.4f}")

    nucleosides = {
        'Adenosine': 'C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)CO)O)O',
        'Guanosine': 'NC1=NC(=O)C2=NCN=C2N1C3OC(CO)C(O)C3O',
        'Cytidine': 'NC1=NC(=O)N=CN1C2OC(CO)C(O)C2O',
        'Uridine': 'O=C1NC=CN=C1C2OC(CO)C(O)C2O'
    }

    print("\nCalculating Balaban J index for the four nucleosides...")
    closest_molecule_name = ""
    smallest_diff = float('inf')

    for name, smi in nucleosides.items():
        mol = Chem.MolFromSmiles(smi)
        balaban_j = GraphDescriptors.BalabanJ(mol)
        diff = abs(balaban_j - ref_balaban_j)
        print(f"  - Balaban J for {name}: {balaban_j:.4f} (Difference: {diff:.4f})")
        if diff < smallest_diff:
            smallest_diff = diff
            closest_molecule_name = name
    
    target_substance_name = closest_molecule_name
    target_substance_smi = nucleosides[target_substance_name]
    print(f"\nThe substance with the most similar Balaban J index is '{target_substance_name}'. This is our target substance.")


    # --- Step 3: Calculate the required indices for the target substance ---
    print(f"\n--- Step 3: Calculating Indices for {target_substance_name} ---")

    target_mol = Chem.MolFromSmiles(target_substance_smi)
    
    # Calculate Zagreb(1) index
    zagreb_m1 = GraphDescriptors.Zagreb1(target_mol)
    print(f"Zagreb(1) index for {target_substance_name}: {zagreb_m1}")

    # Calculate Hosoya Z index (H-included)
    print("Preparing to calculate Hosoya Z index (this may take a moment)...")
    mol_with_H = Chem.AddHs(target_mol)
    
    G = networkx.Graph()
    for atom in mol_with_H.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_with_H.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    
    print(f"Generated H-included graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    
    match_poly = matching_polynomial(G)
    hosoya_z = sum(match_poly.coeffs)
    print(f"Hosoya Z (H-included) index for {target_substance_name}: {int(hosoya_z)}")


    # --- Step 4: Compute the final ratio ---
    print("\n--- Step 4: Computing the Final Ratio ---")
    
    final_ratio = (2 * hosoya_z) / zagreb_m1

    print(f"The final calculation is: (2 * Hosoya Z) / Zagreb(1)")
    print(f"Final Equation: (2 * {int(hosoya_z)}) / {zagreb_m1} = {final_ratio}")
    
    return final_ratio

if __name__ == "__main__":
    final_answer = solve_cheminformatics_problem()
    # The final output requested by the format specification
    print(f"\n<<<{final_answer}>>>")