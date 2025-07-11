import sys
import io

# Set stdout to utf-8 to ensure chemical symbols print correctly
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def solve_wittig_reaction():
    """
    This script determines the product of a Wittig reaction between pivalaldehyde
    and (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane.
    """
    # Step 1: Define reactants
    aldehyde_name = "Pivalaldehyde (2,2-dimethylpropanal)"
    aldehyde_structure = "(CH₃)₃C-CHO"
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_structure = "(o-Cl-C₆H₄)-CH₂-CH=PPh₃"
    byproduct_name = "Triphenylphosphine oxide"
    byproduct_structure = "O=PPh₃"

    print("--- Wittig Reaction Analysis ---")
    print(f"Reactant 1 (Aldehyde): {aldehyde_name}")
    print(f"Reactant 2 (Ylide):      {ylide_name}\n")

    # Step 2: Explain the reaction logic
    print("The Wittig reaction creates an alkene by swapping the carbonyl oxygen of the aldehyde")
    print("with the carbon group of the ylide. A new C=C double bond is formed.\n")

    # Step 3: Identify the combining fragments
    aldehyde_fragment = "(CH₃)₃C-CH="
    ylide_fragment = "=CH-CH₂-(o-Cl-C₆H₄)"

    print("Breaking down the reactants into fragments for combination:")
    print(f"From {aldehyde_name}:  {aldehyde_fragment}")
    print(f"From {ylide_name}: {ylide_fragment}\n")

    # Step 4: Construct the product
    # The ylide is non-stabilized, leading to the Z-isomer as the major product.
    product_structure = aldehyde_fragment + ylide_fragment.lstrip('=')
    product_iupac_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    # SMILES is a standard machine-readable format for molecules.
    # The syntax /C=C\ indicates a Z (cis) double bond.
    product_smiles = "C(C)(C)C/C=C\\Cc1c(Cl)cccc1"

    print("--- Reaction Product Details ---")
    print(f"Combined Structure: {product_structure}")
    print(f"IUPAC Name:         {product_iupac_name}")
    print(f"SMILES String:      {product_smiles}\n")
    
    # Step 5: Display the final reaction equation
    print("--- Final Reaction Equation ---")
    # The term "number" in the prompt's instructions is interpreted as showing the full named equation.
    reactant1 = aldehyde_name
    reactant2 = ylide_name
    final_product = product_iupac_name
    byproduct = byproduct_name
    
    print(f"1. By Name:")
    print(f"  {reactant1}  +  {reactant2}  -->  {final_product}  +  {byproduct}\n")
    
    print(f"2. By Structure:")
    print(f"  {aldehyde_structure}  +  {ylide_structure}  -->  {product_structure}  +  {byproduct_structure}")

# Execute the function to print the solution
solve_wittig_reaction()
