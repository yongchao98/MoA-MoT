import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors

def solve_chemistry_riddle():
    """
    This function solves the chemistry riddle by identifying molecules,
    calculating their graph-theoretical properties, and finding the final product.
    Note: This script requires the rdkit-pypi and mordred packages.
    You can install them using: pip install rdkit-pypi mordred
    """

    # Step 1 & 2: Define molecules via SMILES and create rdkit molecule objects.
    # The deciphering is based on interpreting the riddle's clues as pointing
    # to fundamental biological and aromatic compounds.
    molecules_data = {
        "Guanine (Y1)": "N1C=NC2=C1C(=O)N=C(N)N2",
        "Adenine (Y2)": "NC1=C2N=CNC2=N1",
        "Cytosine (Y3)": "NC1=NC(=O)C=CN1",
        "Thymine (Y4)": "CC1=CN=C(=O)NC1=O",
        "Uracil (Y5)": "C1=CN=C(=O)NC1=O",
        "Phenol (Y6)": "c1ccc(O)cc1",
        "Toluene (Y7)": "Cc1ccccc1",
        "Benzene (Y8)": "c1ccccc1",
        "Naphthalene (Y9)": "c1cccc2ccccc12",
        "Anthracene (Y10)": "c1ccc2cc3ccccc3cc2c1",
    }

    mols = {}
    for name, smi in molecules_data.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols[name] = Chem.AddHs(mol)

    # Step 3: Calculate the Mass-Weighted Barysz Graph Energy for all molecules.
    # The corresponding descriptor in mordred is 'SpAbs_B(m)'.
    barysz_descriptor = descriptors.BaryszMatrix.Barysz(prop='mass')
    calc_barysz = Calculator([barysz_descriptor])
    
    barysz_energies = {name: float(calc_barysz(mol)[0]) for name, mol in mols.items()}

    # Step 4: Find the molecule with the lowest Barysz energy.
    min_energy = float('inf')
    min_energy_molecule_name = None

    for name, energy in barysz_energies.items():
        if energy < min_energy:
            min_energy = energy
            min_energy_molecule_name = name
    
    print(f"Identified molecule (Y) with the lowest energy: {min_energy_molecule_name}")

    # Step 5: For the identified molecule, calculate Mass-Weighted Moran's I for lags 1-8.
    target_mol = mols[min_energy_molecule_name]
    
    moran_descriptors = [descriptors.Moran.Moran(prop='mass', lag=i) for i in range(1, 9)]
    calc_moran = Calculator(moran_descriptors)
    moran_results = calc_moran(target_mol)
    
    # Filter for valid, finite results.
    valid_moran_values = [
        float(val) for val in moran_results 
        if not isinstance(val, (str, Exception)) and np.isfinite(float(val))
    ]
    
    min_moran_I = min(valid_moran_values)
    max_moran_I = max(valid_moran_values)

    # Step 6: Calculate the final product and print the full equation.
    final_product = min_energy * min_moran_I * max_moran_I

    print(f"Lowest Mass-Weighted Barysz Graph Energy: {min_energy}")
    print(f"Minimum Mass-Weighted Moran's I: {min_moran_I}")
    print(f"Maximum Mass-Weighted Moran's I: {max_moran_I}")
    
    print(f"\nFinal Equation:")
    print(f"{min_energy} * {min_moran_I} * {max_moran_I} = {final_product}")

    # Return final answer in the required format
    print(f"\n<<<{final_product}>>>")

if __name__ == '__main__':
    solve_chemistry_riddle()