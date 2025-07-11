# The user might need to install rdkit: pip install rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    print("This script requires the rdkit library.")
    print("Please install it using: pip install rdkit")
    exit()

def solve_wittig_reaction():
    """
    This function determines the product of a specific Wittig reaction,
    and prints the balanced equation with an explanation of the molecular formulas.
    """
    # 1. Define the reactants using SMILES (Simplified Molecular Input Line Entry System)
    # Pivalaldehyde: (CH3)3C-CHO
    pivalaldehyde_smiles = "CC(C)(C)C=O"
    # Ylide: (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane, Ph3P=CH-CH2-(o-Cl-Ph)
    wittig_reagent_smiles = "c1(Cl)ccccc1CC=[P](c2ccccc2)(c3ccccc3)c4ccccc4"

    aldehyde = Chem.MolFromSmiles(pivalaldehyde_smiles)
    ylide = Chem.MolFromSmiles(wittig_reagent_smiles)

    if not aldehyde or not ylide:
        print("Error: Could not parse one or both of the reactant SMILES strings.")
        return

    # 2. Define the Wittig reaction using reaction SMARTS
    # This maps the atoms from reactants to products: C=O + C=P -> C=C + O=P
    # We only care about the main organic product for this reaction.
    rxn_smarts = "[#6:1]=[O:2].[#6:3]=[P]>>[#6:1]=[#6:3]"
    rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)

    # 3. Run the reaction
    reactants = (aldehyde, ylide)
    products = rxn.RunReactants(reactants)

    if not products:
        print("Reaction simulation failed.")
        return

    # 4. Get the product molecule and calculate formulas for the equation
    product_mol = products[0][0]
    Chem.SanitizeMol(product_mol)

    # Calculate formulas for all species in the balanced equation
    aldehyde_formula = "C5H10O"
    ylide_formula = "C26H22PCl"
    alkene_formula = rdMolDescriptors.CalcMolFormula(product_mol)
    phosphine_oxide_formula = "C18H15PO"
    
    # Manually determined IUPAC name for the main product
    product_iupac_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"

    # 5. Print the results as a chemical equation and explain the "numbers"
    print("The main organic product of the reaction is:")
    print(f"IUPAC Name: {product_iupac_name}")

    print("\n---")

    print("The balanced chemical equation for the reaction is:")
    print(f"{aldehyde_formula} + {ylide_formula} -> {alkene_formula} + {phosphine_oxide_formula}")
    
    print("\nHere is a breakdown of each molecule's atomic composition (the numbers in the equation):")
    
    print("\nReactants:")
    print(f"1. Pivalaldehyde ({aldehyde_formula})")
    print("   - C (Carbon): 5 atoms")
    print("   - H (Hydrogen): 10 atoms")
    print("   - O (Oxygen): 1 atom")
    
    print(f"\n2. (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane ({ylide_formula})")
    print("   - C (Carbon): 26 atoms")
    print("   - H (Hydrogen): 22 atoms")
    print("   - P (Phosphorus): 1 atom")
    print("   - Cl (Chlorine): 1 atom")

    print("\nProducts:")
    print(f"1. {product_iupac_name} ({alkene_formula})")
    print("   - C (Carbon): 13 atoms")
    print("   - H (Hydrogen): 17 atoms")
    print("   - Cl (Chlorine): 1 atom")
    
    print(f"\n2. Triphenylphosphine oxide ({phosphine_oxide_formula})")
    print("   - C (Carbon): 18 atoms")
    print("   - H (Hydrogen): 15 atoms")
    print("   - P (Phosphorus): 1 atom")
    print("   - O (Oxygen): 1 atom")


if __name__ == "__main__":
    solve_wittig_reaction()
<<<1-(2-chlorophenyl)-4,4-dimethylpent-2-ene>>>