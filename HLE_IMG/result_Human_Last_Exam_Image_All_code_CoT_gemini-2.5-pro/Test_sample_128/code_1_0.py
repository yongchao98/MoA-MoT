# First, ensure rdkit is installed: pip install rdkit-pypi
from rdkit import Chem
from rdkit.Chem import Descriptors

def identify_product_A():
    """
    Analyzes the given chemical reaction and identifies the final product, Compound A.
    The code then prints the properties of this compound.
    """
    # Define reaction parameters based on the user's query
    reactant_1_name = "2-aminopyridine"
    reactant_2_name = "o-phthalaldehyde"
    reactant_3_name = "TMSCN"
    temperature_c = 28
    duration_h = 4
    
    print("--- Reaction Analysis ---")
    print(f"The reaction involves {reactant_1_name}, {reactant_2_name}, and {reactant_3_name}.")
    print(f"It is conducted at {temperature_c}°C for {duration_h} hours.")
    print("This is a three-component reaction leading to a stable heterocyclic product.")
    
    # Define the final product based on chemical principles
    # The reaction is a Strecker synthesis followed by intramolecular cyclization and dehydration.
    product_A_name = "2-(pyridin-2-yl)isoindole-1-carbonitrile"
    product_A_smiles = "N#Cc1c2ccccc2n1c3ncccc3" # SMILES representation of the product

    # Use RDKit to process the chemical information
    mol = Chem.MolFromSmiles(product_A_smiles)
    
    if mol:
        molecular_formula = Descriptors.MolFormula(mol)
        molecular_weight = Descriptors.MolWt(mol)
        
        print("\n--- Product A Details ---")
        print(f"The structure of compound A is: {product_A_name}")
        print(f"Molecular Formula: {molecular_formula}")
        print(f"Molecular Weight: {molecular_weight:.2f}")
        print(f"SMILES String: {product_A_smiles}")

        # As requested, outputting the numbers from the problem description
        print("\n--- Numbers from Reaction Description ---")
        print(f"From '2-aminopyridine': 2")
        print(f"From '28 C': {temperature_c}")
        print(f"From '4 h': {duration_h}")
    else:
        print("Error: Could not generate molecule from SMILES string.")

# Run the analysis
identify_product_A()