import re
try:
    import pubchempy as pcp
except ImportError:
    print("PubChemPy library not found. Please install it using: pip install pubchempy")
    pcp = None

# --- Chemical Explanation ---
print("### Chemical Reaction Analysis ###")
print("The reaction proceeds in two steps:\n")
print("Step 1: Directed ortho-Metalation (DoM)")
print("The starting material, N,N-diethyl-3-dimethylaminobenzamide, has two directing groups:")
print("  - The amide group (-CONEt2) at position 1, which directs lithiation to position 2.")
print("  - The dimethylamino group (-NMe2) at position 3, which also directs lithiation to position 2.")
print("The strong base (sec-BuLi) removes the proton from position 2, forming an aryllithium intermediate.\n")

print("Step 2: Electrophilic Quench with Methyl Iodide")
print("The aryllithium intermediate is a potent nucleophile that reacts with methyl iodide (CH3I).")
print("The lithium atom at position 2 is replaced by a methyl group.\n")

# --- Product Identification ---
reactant_name = "N,N-diethyl-3-dimethylaminobenzamide"
product_name_predicted = "N,N-diethyl-2-methyl-3-(dimethylamino)benzamide"
reagents = "1) sec-BuLi/TMEDA, 2) CH3I"
final_product_name = ""

# The SMILES (Simplified Molecular-Input Line-Entry System) string for the predicted product
product_smiles = "CN(C)C1=C(C)C(=CC=C1)C(=O)N(CC)CC"

if pcp:
    try:
        # Use PubChem to find the canonical IUPAC name from the SMILES string
        compounds = pcp.get_compounds(product_smiles, 'smiles')
        if compounds:
            final_product_name = compounds[0].iupac_name
        else:
            final_product_name = product_name_predicted
    except Exception:
        # Fallback in case of network issues
        final_product_name = product_name_predicted
else:
    final_product_name = product_name_predicted

print("### Final Result ###")
print(f"The final compound obtained is: {final_product_name}")
print("\nThe full reaction can be written as:")
print(f"{reactant_name}  ---[{reagents}]--->  {final_product_name}\n")

# --- Extracting Numbers from the Equation ---
reactant_numbers = re.findall(r'\d+', reactant_name)
product_numbers = re.findall(r'\d+', final_product_name)
step_numbers = ['1', '2'] # From the reagents string

print("As requested, here are the numbers from the reaction equation:")
print(f"Number in the reactant name '{reactant_name}': {', '.join(reactant_numbers)}")
print(f"Numbers for the reaction steps: {', '.join(step_numbers)}")
print(f"Numbers in the product name '{final_product_name}': {', '.join(product_numbers)}")

<<<N,N-diethyl-2-methyl-3-(dimethylamino)benzamide>>>