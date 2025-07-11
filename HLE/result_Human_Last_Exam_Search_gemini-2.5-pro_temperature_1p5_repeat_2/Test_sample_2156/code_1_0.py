import pandas as pd
from rdkit import Chem
from mordred import Calculator, descriptors

def solve():
    """
    Solves the problem by calculating molecular descriptors for formaldehyde homologs.
    """
    aldehydes = {
        'Acetaldehyde': 'CC=O',
        'Propionaldehyde': 'CCC=O',
        'Butyraldehyde': 'CCCC=O',
        'Pentanal': 'CCCCC=O',
        'Hexanal': 'CCCCCC=O',
        'Heptanal': 'CCCCCCC=O',
        'Octanal': 'CCCCCCCC=O',
        'Nonanal': 'CCCCCCCCC=O',
        'Decanal': 'CCCCCCCCCC=O',
    }

    # Initialize Mordred calculator
    calc = Calculator(descriptors, ignore_3D=True)

    # We are interested in GATS (Geary autocorrelation) weighted by Sanderson electronegativity ('se')
    # and average path chi indices.
    calc.register([
        descriptors.Autocorrelation.GATSse,
        descriptors.Chi.Path,
        descriptors.Chi.AveragePath,
    ])


    results = []

    for name, smiles in aldehydes.items():
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        # Calculate all descriptors
        try:
            all_descriptors = calc(mol)
        except Exception:
            continue

        # Extract GATSse values for lags 1 to 8
        gats_values = {}
        for i in range(1, 9):
            try:
                # Descriptor names are like GATS1se, GATS2se etc.
                gats_values[i] = all_descriptors[f'GATS{i}se']
            except KeyError:
                # If a descriptor is not calculated (e.g., due to molecule size), skip it
                continue

        if not gats_values:
            continue

        # Find the maximum GATS value and the corresponding lag (imax)
        imax = max(gats_values, key=gats_values.get)
        max_gats = gats_values[imax]

        # Check if the max GATS value is between 2 and 3
        if 2 <= max_gats <= 3:
            # Calculate average simple and valence path chi indices
            # ASP (Average Simple Path) uses 'd' property
            # AVP (Average Valence Path) uses 'dv' property
            avg_simple_chi_values = []
            avg_valence_chi_values = []
            for i in range(11): # Orders 0 to 10
                 # Mordred descriptor names are like AXp-0d, AXp-1d etc for simple path
                 # and AXp-0dv, AXp-1dv etc for valence path.
                try:
                    simple_val = all_descriptors[f'AXp-{i}d']
                    valence_val = all_descriptors[f'AXp-{i}dv']
                    if simple_val is not None and not pd.isna(simple_val) and simple_val > 0 :
                        avg_simple_chi_values.append(simple_val)
                    if valence_val is not None and not pd.isna(valence_val) and valence_val > 0:
                        avg_valence_chi_values.append(valence_val)
                except KeyError:
                    # Not all orders might be calculated for all molecules
                    continue

            if not avg_simple_chi_values or not avg_valence_chi_values:
                continue

            # Calculate the average of the non-zero values for each index
            avg_simple_chi = sum(avg_simple_chi_values) / len(avg_simple_chi_values)
            avg_valence_chi = sum(avg_valence_chi_values) / len(avg_valence_chi_values)

            # Calculate the difference
            chi_difference = avg_valence_chi - avg_simple_chi

            # Calculate the final product
            product = imax * chi_difference

            results.append({
                'homolog': name,
                'imax': imax,
                'max_gats': max_gats,
                'chi_difference': chi_difference,
                'product': product,
            })

    if not results:
        # No homologs were found in the specified GATS range.
        # This is a possible outcome depending on the data.
        return None

    # Find the minimum product among the filtered homologs
    min_result = min(results, key=lambda x: x['product'])

    return min_result['product']

# In a functional environment, this would be called to get the answer.
# final_answer = solve()
# print(f"<<<{final_answer}>>>")