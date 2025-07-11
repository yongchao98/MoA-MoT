# The user needs to have rdkit-pypi, mordred, and pandas installed.
# If not installed, they can be installed using pip:
# pip install rdkit-pypi mordred pandas numpy

from rdkit import Chem
from mordred import Calculator
from mordred.GearyAuto import GearyAuto
from mordred.PathCount import AvcPath, Sxp
import pandas as pd
import numpy as np
import sys

def find_minimum_product_of_homologs():
    """
    This function executes the plan to find the minimum product based on
    molecular descriptors for formaldehyde's homologs.
    """
    # Step 1: Define homologs of formaldehyde (aldehydes C1-C12).
    # Using a list of names and their corresponding SMILES strings.
    aldehydes = {
        "formaldehyde": "C=O", "acetaldehyde": "CC=O", "propanal": "CCC=O",
        "butanal": "CCCC=O", "pentanal": "CCCCC=O", "hexanal": "CCCCCC=O",
        "heptanal": "CCCCCCC=O", "octanal": "CCCCCCCC=O", "nonanal": "CCCCCCCCC=O",
        "decanal": "CCCCCCCCCC=O", "undecanal": "CCCCCCCCCCC=O", "dodecanal": "CCCCCCCCCCCC=O"
    }

    # Create RDKit molecule objects, adding hydrogens as they are important for some calculations.
    try:
        mols = {name: Chem.AddHs(Chem.MolFromSmiles(smi)) for name, smi in aldehydes.items()}
    except Exception as e:
        print(f"Error creating molecules with RDKit: {e}")
        return

    # Step 2: Set up the mordred descriptor calculator.
    calc = Calculator(ignore_3D=True)
    
    # Register Geary autocorrelations weighted by Sanderson electronegativity ('pe') for lags 1-8.
    for i in range(1, 9):
        calc.register(GearyAuto(i, 'pe'))
    
    # Register the average chi path indices.
    calc.register(AvcPath)
    calc.register(Sxp)

    # Calculate descriptors for all molecules and store results in a pandas DataFrame.
    # The 'quiet=True' argument suppresses the calculation progress bar from appearing in the output.
    df = calc.pandas(mols.values(), quiet=True)
    df.index = mols.keys()

    # Handle potential calculation failures by ensuring all descriptor columns exist.
    gats_cols = [f'GATS{i}pe' for i in range(1, 9)]
    for col in gats_cols + ['AvcPath', 'Sxp']:
        if col not in df.columns:
            df[col] = np.nan # Fill missing columns with Not-a-Number
    df.dropna(inplace=True) # Drop rows where any descriptor failed

    # Step 3: Find the maximum GATS value and the corresponding lag (i_max) for each homolog.
    df_gats = df[gats_cols]
    df['max_GATS'] = df_gats.max(axis=1)
    df['i_max_col'] = df_gats.idxmax(axis=1)
    
    # Extract the lag 'i' from the column name (e.g., from 'GATS3pe', extract 3).
    df['i_max'] = df['i_max_col'].str.extract(r'GATS(\d+)pe').astype(int)

    # Step 4: Filter the homologs, keeping only those where 2 < max_GATS < 3.
    filtered_homologs = df[(df['max_GATS'] > 2) & (df['max_GATS'] < 3)].copy()

    # Check if any homologs satisfied the condition.
    if filtered_homologs.empty:
        print("No homologs were found that satisfy the condition on Geary autocorrelation (2 < max_GATS < 3).")
        print("Therefore, the minimum product cannot be determined.")
        return

    # Step 5: For the filtered homologs, calculate the target product.
    filtered_homologs['diff_chi'] = filtered_homologs['AvcPath'] - filtered_homologs['Sxp']
    filtered_homologs['product'] = filtered_homologs['i_max'] * filtered_homologs['diff_chi']

    # Step 6: Identify the homolog with the minimum product.
    min_product_row = filtered_homologs.loc[filtered_homologs['product'].idxmin()]

    # Step 7: Print the detailed results and the final answer.
    min_product_val = min_product_row['product']
    i_max_val = min_product_row['i_max']
    avcpath_val = min_product_row['AvcPath']
    sxp_val = min_product_row['Sxp']
    homolog_name = min_product_row.name

    print(f"The minimum product was found for the homolog: {homolog_name}")
    print("The final equation is the product of i_max and the difference between average chi indices:")
    print(f"product = {i_max_val} * ({avcpath_val:.5f} - {sxp_val:.5f})")
    
    # The final answer in the required format.
    print(f"\n<<<{min_product_val}>>>")


if __name__ == '__main__':
    # Checking for required packages before running the main function.
    try:
        import rdkit
        import mordred
        import pandas
        import numpy
    except ImportError:
        print("A required library (rdkit, mordred, pandas, or numpy) is not installed.")
        print("Please install them using: pip install rdkit-pypi mordred pandas numpy")
        sys.exit(1)
    
    find_minimum_product_of_homologs()