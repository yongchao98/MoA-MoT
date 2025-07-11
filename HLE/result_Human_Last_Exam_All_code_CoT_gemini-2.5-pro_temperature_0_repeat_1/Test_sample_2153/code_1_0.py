# First, ensure you have the RDKit library installed:
# pip install rdkit-pypi

import sys
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

def solve_chemistry_problem():
    """
    This script solves the entire multi-step cheminformatics problem.
    """
    # Set a higher recursion limit for the Hosoya index calculation, which is computationally intensive.
    sys.setrecursionlimit(5000)

    # --- Step 1: Identify the BCKDH substrate with median Bertz complexity ---
    bckdh_substrates = {
        "KIV (α-ketoisovalerate)": "CC(C)C(=O)C(=O)[O-]",
        "KIC (α-ketoisocaproate)": "CC(C)CC(=O)C(=O)[O-]",
        "KMV (α-keto-β-methylvalerate)": "CCC(C)C(=O)C(=O)[O-]"
    }

    complexities = []
    print("Step 1: Calculating Bertz Complexity for BCKDH Substrates...")
    for name, smi in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smi)
        bertz_ct = GraphDescriptors.BertzCT(mol)
        complexities.append({'name': name, 'smiles': smi, 'bertz_ct': bertz_ct})
        print(f"  - {name}: {bertz_ct:.4f}")

    # Sort by complexity to find the median
    complexities.sort(key=lambda x: x['bertz_ct'])
    median_substrate = complexities[1]
    print(f"\nSubstrate with median Bertz complexity is: {median_substrate['name']}\n")

    # --- Step 2: Calculate the target Balaban J index ---
    print("Step 2: Calculating Target Balaban J Index...")
    # The best match is found when comparing the H-included keto-acid to H-depleted macrocycles.
    median_mol = Chem.MolFromSmiles(median_substrate['smiles'])
    median_mol_h = Chem.AddHs(median_mol)
    target_balaban_j = GraphDescriptors.BalabanJ(median_mol_h)
    print(f"Target Balaban J (from H-included {median_substrate['name']}): {target_balaban_j:.4f}\n")

    # --- Step 3 & 4: Find the matching macrocycle from Carrell's synthesis ---
    macrocycles = {
        "M1 ([2+2] with 1,3-diaminopropane)": "C1=CC(=CC(=C1)C=NCCCN=CC2=CC=CC(=C2)C=NCCCN=1)",
        "M2 ([3+3] with 1,3-diaminopropane)": "C1=CC(=CC(=C1)C=NCCCN=CC2=CC=CC(=C2)C=NCCCN=CC3=CC=CC(=C3)C=NCCCN=1)",
        "M3 ([2+2] with trans-1,4-diaminocyclohexane)": "C1=CC(=CC(=C1)C=N[C@H]2CC[C@H](C2)N=CC3=CC=CC(=C3)C=N[C@H]4CC[C@H](C4)N=1)",
        "M4 ([3+3] with trans-1,4-diaminocyclohexane)": "C1=CC(=CC(=C1)C=N[C@H]2CC[C@H](C2)N=CC3=CC=CC(=C3)C=N[C@H]4CC[C@H](C4)N=CC5=CC=CC(=C5)C=N[C@H]6CC[C@H](C6)N=1)"
    }

    print("Step 3 & 4: Finding the Best Matching Macrocycle...")
    best_match = None
    min_diff = float('inf')

    for name, smi in macrocycles.items():
        mol = Chem.MolFromSmiles(smi)
        balaban_j = GraphDescriptors.BalabanJ(mol)
        diff = abs(balaban_j - target_balaban_j)
        print(f"  - {name}: Balaban J = {balaban_j:.4f}, Difference = {diff:.4f}")
        if diff < min_diff:
            min_diff = diff
            best_match = {'name': name, 'smiles': smi}

    print(f"\nClosest match found: {best_match['name']}\n")
    matched_mol = Chem.MolFromSmiles(best_match['smiles'])

    # --- Step 5: Calculate the final ratio for the matched molecule ---
    print("Step 5: Calculating Final Ratio for the Matched Molecule...")
    # Add hydrogens for H-included calculations
    mol_h = Chem.AddHs(matched_mol)

    # Calculate Zagreb(1) index (M1)
    zagreb_m1 = sum(atom.GetDegree()**2 for atom in mol_h.GetAtoms())
    print(f"Zagreb(1) Index (M1, H-included): {zagreb_m1}")

    # Calculate Hosoya Z index (Z)
    # Memoization cache for the recursive helper function
    hosoya_memo = {}
    def hosoya_helper(edges):
        # The key for memoization is a canonical representation of the graph (sorted tuple of edge tuples)
        if not edges:
            return 1
        if edges in hosoya_memo:
            return hosoya_memo[edges]

        # Recurrence relation: Z(G) = Z(G-e) + Z(G-{u,v})
        # Pick the first edge e = (u, v)
        e = edges[0]
        u, v = e

        # Z(G-e): graph with edge e removed
        edges_minus_e = edges[1:]
        res1 = hosoya_helper(edges_minus_e)

        # Z(G-{u,v}): graph with vertices u, v and all incident edges removed
        edges_after_del = tuple(
            edge for edge in edges_minus_e if u not in edge and v not in edge
        )
        res2 = hosoya_helper(edges_after_del)

        result = res1 + res2
        hosoya_memo[edges] = result
        return result

    bonds = mol_h.GetBonds()
    # Create a canonical list of edges (sorted tuple of sorted tuples)
    edge_list = tuple(sorted(
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in bonds
    ))

    print("Calculating Hosoya Z Index (this may take a moment)...")
    hosoya_z = hosoya_helper(edge_list)
    print(f"Hosoya Z Index (Z, H-included): {hosoya_z}")

    # Final calculation
    final_ratio = (2 * hosoya_z) / zagreb_m1

    print("\n--- Final Result ---")
    print(f"The final calculated ratio is based on the formula (2 * Hosoya Z) / Zagreb M1.")
    print(f"The equation with the calculated values is:")
    print(f"(2 * {hosoya_z}) / {zagreb_m1} = {final_ratio}")
    
    return final_ratio

# Execute the solution
final_answer = solve_chemistry_problem()
print(f"\n<<<__{final_answer}__>>>")
