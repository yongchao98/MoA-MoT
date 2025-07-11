import sys
import warnings
from rdkit import Chem
# The mordred library can be installed via pip: pip install mordred
try:
    from mordred.GearyAuto import GearyAuto
    from mordred.Chi import AAvgValencePath, AAvgSimplePath
except ImportError:
    print("Mordred library not found. Please install it using: pip install mordred")
    sys.exit(1)

# Suppress warnings from RDKit for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

def solve_task():
    """
    Solves the described problem by finding the minimum product of i_max and a difference of chi indices
    for a specific subset of formaldehyde's homologs.
    """
    aldehydes = {
        # C1
        'Formaldehyde': 'C=O',
        # C2
        'Acetaldehyde': 'CC=O',
        # C3
        'Propanal': 'CCC=O',
        # C4
        'Butanal': 'CCCC=O',
        'Isobutanal': 'CC(C)C=O',
        # C5
        'Pentanal': 'CCCCC=O',
        '2-Methylbutanal': 'CCC(C)C=O',
        '3-Methylbutanal': 'CC(C)CC=O',
        'Pivaldehyde': 'C(C)(C)C=O',
        # C6
        'Hexanal': 'CCCCCC=O',
        '2-Methylpentanal': 'CCCC(C)C=O',
        '3-Methylpentanal': 'CCC(C)CC=O',
        '4-Methylpentanal': 'CC(C)CCC=O',
        '2,2-Dimethylbutanal': 'CCC(C)(C)C=O',
        '2,3-Dimethylbutanal': 'CC(C)C(C)C=O',
        '3,3-Dimethylbutanal': 'CC(C)(C)CC=O'
    }

    found_homologs_data = []

    for name, smiles in aldehydes.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)

        max_geary = -float('inf')
        i_max = -1

        # Calculate Geary autocorrelation for lags 1 through 8
        for i in range(1, 9):
            try:
                geary_desc = GearyAuto(lag=i, prop='electronegativity_sanderson')
                geary_val = geary_desc(mol)
                if not isinstance(geary_val, (float, int)):
                    continue
                if geary_val > max_geary:
                    max_geary = geary_val
                    i_max = i
            except Exception:
                continue

        # Check if the homolog meets the criterion
        if 2 < max_geary < 3:
            try:
                avg_val_path = float(AAvgValencePath()(mol))
                avg_sim_path = float(AAvgSimplePath()(mol))
                
                diff = avg_val_path - avg_sim_path
                product = i_max * diff
                
                found_homologs_data.append({
                    'name': name,
                    'i_max': i_max,
                    'diff': diff,
                    'product': product
                })
            except Exception:
                continue

    if not found_homologs_data:
        print("No homologs were found that match the specified Geary autocorrelation criteria.")
        return

    # Find the homolog with the minimum product
    min_homolog = min(found_homologs_data, key=lambda x: x['product'])

    # Print the final result in the format of an equation
    print(f"{min_homolog['i_max']} * ({float(AAvgValencePath()(Chem.AddHs(Chem.MolFromSmiles(aldehydes[min_homolog['name']]))))} - {float(AAvgSimplePath()(Chem.AddHs(Chem.MolFromSmiles(aldehydes[min_homolog['name']]))))}) = {min_homolog['product']}")


if __name__ == "__main__":
    solve_task()