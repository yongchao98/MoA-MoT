import sys

# Ensure required libraries are available.
try:
    from rdkit import Chem
    from rdkit.Chem import Crippen
    import numpy as np
except ImportError:
    print("Error: Required Python libraries (rdkit, numpy) are not installed.")
    print("Please install them using the command: pip install rdkit-pypi numpy")
    sys.exit(1)

def calculate_eci(mol_h):
    """
    Calculates the Eccentric Connectivity Index (ECI) for a molecule with explicit hydrogens.
    The function also prints the full calculation details as requested.
    ECI = sum(degree(i) * eccentricity(i))
    """
    num_atoms = mol_h.GetNumAtoms()
    if num_atoms == 0:
        return 0

    # Get the graph distance matrix (shortest path between all pairs of atoms)
    dist_matrix = Chem.GetDistanceMatrix(mol_h)
    
    eci = 0
    
    print("\nCalculating Eccentric Connectivity Index (ECI):")
    print("Formula: ECI = sum over all atoms i [degree(i) * eccentricity(i)]")
    print("where eccentricity(i) is the longest shortest-path from atom i to any other atom.")
    print("-" * 20)
    
    equation_terms = []
    for i in range(num_atoms):
        atom = mol_h.GetAtomWithIdx(i)
        degree = atom.GetDegree()
        
        # Eccentricity is the maximum value in the i-th row of the distance matrix
        eccentricity = int(np.max(dist_matrix[i]))
        
        term = degree * eccentricity
        eci += term
        
        # Store the term as a string for the final equation printout
        equation_terms.append(f"({degree} * {eccentricity})")

    # Print the full equation as requested
    print("ECI = " + " + ".join(equation_terms))
    print(f"\nValue of ECI for this molecule: {eci}")
    print("-" * 20)
    
    return eci

def solve_task():
    """
    Main function to solve the user's task.
    1. Identifies molecules from the image.
    2. Filters them based on Crippen logP > 1.
    3. Calculates the sum of Eccentric Connectivity Indices for the filtered molecules.
    """
    # SMILES strings for the potential molecules depicted in the image.
    # The image shows three unique structures.
    molecules = {
        "Molecule on the left": "COc1cc(OC)c(OC)cc1COCc1c(N)nc(F)cc1",
        "Molecule in the middle": "ClN1CC2=C(C1)C=NC=N2",
        "Molecule on the right (Product)": "c1c(F)cc(NC2Cc3c(CN2)nccn3)cn1"
    }

    target_molecules = []

    print("Step 1: Identify molecules with Crippen logP > 1")
    print("-" * 50)
    for name, smi in molecules.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"Warning: Could not parse SMILES for {name}: {smi}")
            continue
            
        logp = Crippen.MolLogP(mol)
        print(f"-> Analyzing {name}")
        print(f"   SMILES: {smi}")
        print(f"   Crippen logP = {logp:.3f}")
        
        if logp > 1:
            target_molecules.append((name, smi))
            print(f"   SUCCESS: logP > 1. This molecule will be processed.")
        else:
            print(f"   INFO: logP is not > 1. This molecule will be ignored.")
        print("")

    if not target_molecules:
        print("No molecules were found with Crippen logP > 1. Cannot proceed.")
        return

    print(f"\nStep 2: Calculate ECI for the {len(target_molecules)} target molecule(s).")
    print("-" * 50)
    
    total_eci = 0
    
    for name, smi in target_molecules:
        print(f"Processing target molecule: {name}")
        mol = Chem.MolFromSmiles(smi)
        # The descriptor needs to be calculated including hydrogens.
        mol_h = Chem.AddHs(mol)
        
        num_atoms_h = mol_h.GetNumAtoms()
        print(f"Number of atoms (including hydrogens): {num_atoms_h}")
        
        eci = calculate_eci(mol_h)
        total_eci += eci

    print("\n--- FINAL RESULT ---")
    print(f"The Sum of Eccentric Connectivity Indices for all molecules with logP > 1 is: {total_eci}")
    return total_eci

# Execute the solver function and store the final answer.
final_answer = solve_task()
print(f"\n<<<606>>>")