import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def solve_chemistry_puzzle():
    """
    Solves the multi-step chemistry puzzle by identifying molecules, calculating their
    properties, and determining the final result.
    """

    # Step 1: Define the list of plausible explosive molecules based on the problem's context.
    # The SMILES strings represent the molecular structure.
    molecules_data = [
        ("Y1_TNT", "CC1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]"),
        ("Y2_Nitroglycerin", "C(C(CO[N+](=O)[O-])O[N+](=O)[O-])O[N+](=O)[O-]"),
        ("Y3_RDX", "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]"),
        ("Y4_PETN", "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])O[N+](=O)[O-]"),
        ("Y5_Picric_Acid", "C1=C(C(=C(C=C1[N+](=O)[O-])O)[N+](=O)[O-])[N+](=O)[O-]"),
        # Using a simplified monomer representation for Guncotton (Cellulose Trinitrate)
        ("Y6_Guncotton", "C(C1C(C(C(C1O)O[N+](=O)[O-])O[N+](=O)[O-])O[N+](=O)[O-])O"),
        # DDNP (Diazodinitrophenol)
        ("Y7_DDNP", "c1c([N+](=O)[O-])cc(c(O)c1[N+](=O)[O-])C=[N+]=[N-]"),
        ("Y8_Mercury_Fulminate", "[O-][N+]#C[Hg]C#[N+][O-]"),
        # Lead Azide is ionic, we'll model one formula unit for graph calculations.
        # This is a simplified linear model for calculation purposes. N=N=N-[Pb]-N=N=N
        ("Y9_Lead_Azide", "N=N=N[Pb]N=N=N"),
        # TATP (Acetone Peroxide Trimer)
        ("Y10_TATP", "CC1(C)OOC(C)(C)OOC(C)(C)OO1")
    ]

    results = []

    for name, smiles in molecules_data:
        # Step 2 & 3: Perform calculations for each molecule
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
        
        mol = Chem.AddHs(mol)
        
        num_atoms = mol.GetNumAtoms()
        adj_matrix = Chem.GetAdjacencyMatrix(mol)
        
        # Get atomic masses
        masses = np.array([atom.GetMass() for atom in mol.GetAtoms()])

        # Calculate Mass-Weighted Barysz Matrix
        barysz_matrix = np.zeros((num_atoms, num_atoms))
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                if adj_matrix[i, j] == 1:
                    val = 1.0 / np.sqrt(masses[i] * masses[j])
                    barysz_matrix[i, j] = barysz_matrix[j, i] = val
        
        # Calculate Barysz Energy
        eigenvalues = np.linalg.eigvalsh(barysz_matrix)
        barysz_energy = np.sum(np.abs(eigenvalues))

        # Calculate Local Moran's I values using atomic mass as the property
        mean_mass = np.mean(masses)
        mass_deviations = masses - mean_mass
        
        # Denominator term for full Moran's I, often represented as S0 and M2
        m2 = np.sum(mass_deviations**2)
        if m2 == 0: # Avoid division by zero for single-atom graphs etc.
            min_moran_I = 0
            max_moran_I = 0
        else:
            # Numerator for local moran I
            local_moran_numerators = mass_deviations * (adj_matrix @ mass_deviations)
            local_I = local_moran_numerators / m2 * num_atoms
            min_moran_I = np.min(local_I)
            max_moran_I = np.max(local_I)
        
        results.append({
            "name": name,
            "energy": barysz_energy,
            "min_moran": min_moran_I,
            "max_moran": max_moran_I
        })
        
    # Step 4: Identify the molecule with the lowest Barysz Energy
    if not results:
        print("Could not process any molecules.")
        return

    min_energy_molecule = min(results, key=lambda x: x["energy"])
    
    # Step 5: Final Calculation
    name = min_energy_molecule["name"]
    energy = min_energy_molecule["energy"]
    min_moran = min_energy_molecule["min_moran"]
    max_moran = min_energy_molecule["max_moran"]
    
    product = energy * min_moran * max_moran
    
    print(f"The identified molecule (Y) with the lowest Mass-Weighted Barysz Graph Energy is: {name}")
    print(f"Identified Energy: {energy}")
    print(f"Minimum Mass-Weighted Moran's I: {min_moran}")
    print(f"Maximum Mass-Weighted Moran's I: {max_moran}")
    print("\nThe final product is calculated as:")
    print(f"Product = {energy} * {min_moran} * {max_moran} = {product}")
    print(f"\n<<<{product}>>>")


solve_chemistry_puzzle()