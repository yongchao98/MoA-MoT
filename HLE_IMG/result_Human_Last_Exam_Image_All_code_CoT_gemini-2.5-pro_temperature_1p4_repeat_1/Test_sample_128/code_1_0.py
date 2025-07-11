try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit library not found. Please install it using: pip install rdkit-pypi")
    # Exit gracefully if rdkit is not available
    exit()

# The final product A is 1-cyano-2-(pyridin-2-yl)isoindole
product_name = "1-cyano-2-(pyridin-2-yl)isoindole"

# SMILES (Simplified Molecular Input Line Entry System) string for the product A
product_smiles = "N#Cc1c2ccccc2n1c3ncccc3"

# Create a molecule object from the SMILES string
mol = Chem.MolFromSmiles(product_smiles)

if mol:
    # Calculate properties of the product
    molecular_formula = CalcMolFormula(mol)
    exact_molecular_weight = Descriptors.ExactMolWt(mol)

    # Print the reaction and the identified product
    print("The reaction is a three-component condensation.")
    print("Product A is identified as: {}".format(product_name))
    print("\n--- Chemical Properties of Product A ---")
    print("SMILES String: {}".format(product_smiles))
    print("Molecular Formula: {}".format(molecular_formula))
    print("Exact Molecular Weight: {:.4f}".format(exact_molecular_weight))
    
    # Print the balanced chemical equation with each component
    print("\n--- Balanced Chemical Equation ---")
    print("Reactants:")
    print("  1 molecule of 2-aminopyridine (C5H6N2)")
    print("  1 molecule of o-phthalaldehyde (C8H6O2)")
    print("  1 molecule of Hydrogen Cyanide (HCN) from TMSCN")
    print("Products:")
    print("  1 molecule of Product A ({})".format(molecular_formula))
    print("  2 molecules of Water (H2O)")
    print("\nEquation:")
    print("C5H6N2 + C8H6O2 + HCN -> {} + 2H2O".format(molecular_formula))

else:
    print("Error: Could not parse the SMILES string for the product.")
