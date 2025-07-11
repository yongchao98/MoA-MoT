import sys
try:
    import pubchempy as pcp
except ImportError:
    print("PubChemPy library not found. Please install it to run this script:")
    print("pip install pubchempy")
    sys.exit(1)

# Step 1: Define the reactants and products for the chemical equation.
# The reaction involves 1 equivalent of each reactant.
reactant1 = "Methyl phenyl sulfoxide (C6H5S(O)CH3)"
reactant2 = "Triflic anhydride ((CF3SO2)2O)"
reactant3 = "Trimethylsilyl cyanide ((CH3)3SiCN)"

# The main organic product of this Pummerer reaction is (phenylthio)acetonitrile.
product_main = "2-(Phenylthio)ethanenitrile (C6H5SCH2CN)"
byproduct1 = "Trimethylsilyl triflate ((CF3SO2)OSi(CH3)3)"
byproduct2 = "Triflic acid (CF3SO2OH)"

# Step 2: Print the balanced chemical equation, including the stoichiometric coefficients (the numbers).
print("The balanced chemical equation for the reaction is:\n")
print(
    f"1 {reactant1} + "
    f"1 {reactant2} + "
    f"1 {reactant3}  --->"
)
print(
    f"1 {product_main} + "
    f"1 {byproduct1} + "
    f"1 {byproduct2}"
)
print("\n-------------------------------------------------\n")

# Step 3: Use the product's structure to find its IUPAC name programmatically.
# The SMILES (Simplified Molecular Input Line Entry System) string for the main product is "N#CCSc1ccccc1".
product_smiles = "N#CCSc1ccccc1"

try:
    # Search PubChem for the compound using its SMILES string.
    compounds = pcp.get_compounds(product_smiles, 'smiles')
    
    if compounds:
        # Get the first matching compound.
        compound = compounds[0]
        # Retrieve and print its IUPAC name.
        print(f"The main organic product has the SMILES string: {product_smiles}")
        print(f"The IUPAC name of the product is: {compound.iupac_name}")
    else:
        # Fallback if the compound is not found in PubChem via this script.
        print(f"Could not find a compound with SMILES: {product_smiles} in the PubChem database.")
        print("Based on nomenclature rules, the IUPAC name is '2-(phenylthio)ethanenitrile'.")

except Exception as e:
    print(f"An error occurred while connecting to the PubChem database: {e}")
    print("Based on nomenclature rules, the IUPAC name of C6H5SCH2CN is '2-(phenylthio)ethanenitrile'.")
