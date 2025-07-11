import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors
import networkx as nx

def calculate_eccentric_connectivity_index():
    """
    Identifies the product molecule from the reaction, verifies its Crippen logP,
    and calculates its Eccentric Connectivity Index including hydrogen atoms.
    """
    # Step 1: Define the product molecule using its canonical SMILES.
    # The molecule is N-(5-fluoropyridin-2-yl)-6,7-dihydro-5H-pyrrolo[2,3-d]pyrimidin-4-amine.
    smiles = "C1=C(C=NC(=C1)F)NC2=NC=NC3=C2NCCC3"
    mol = Chem.MolFromSmiles(smiles)

    # Step 2: Verify the Crippen logP > 1 condition.
    logp = Crippen.MolLogP(mol)
    if logp <= 1:
        print(f"The molecule's logP is {logp:.2f}, which does not meet the > 1 criterion.")
        return

    # Step 3: Create the molecular graph including hydrogen atoms.
    mol_h = Chem.AddHs(mol)
    G = nx.Graph()
    for atom in mol_h.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Step 4: Calculate degrees and eccentricities for all atoms.
    num_atoms = mol_h.GetNumAtoms()
    try:
        eccentricities = nx.eccentricity(G)
    except nx.NetworkXError:
        print("Graph is not connected. Cannot compute eccentricity.")
        return
        
    atom_degrees = [val for (node, val) in G.degree()]

    # Step 5: Compute the Eccentric Connectivity Index.
    total_eci = 0
    equation_terms = []
    
    for i in range(num_atoms):
        degree = atom_degrees[i]
        ecc = eccentricities[i]
        term = degree * ecc
        total_eci += term
        equation_terms.append(f"{degree}*{ecc}")

    # Step 6: Print the detailed calculation and the final result.
    print(f"Molecule: {smiles}")
    print(f"Formula: {rdMolDescriptors.CalcMolFormula(mol_h)}")
    print(f"Crippen logP: {logp:.2f}")
    print(f"Total number of atoms (including H): {num_atoms}")
    print("\nCalculating the Eccentric Connectivity Index (XI = sum(degree * eccentricity)):")
    
    full_equation = " + ".join(equation_terms)
    print(f"XI = {full_equation}")
    
    print(f"\nTotal Sum of Eccentric Connectivity Indices = {total_eci}")

calculate_eccentric_connectivity_index()