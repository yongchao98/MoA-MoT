import re
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit library is not installed. Please install it using: pip install rdkit")
    exit()

def get_diels_alder_product():
    """
    This function models the Diels-Alder reaction between butadiene and 
    1,1-dichloro-2,2-difluoroethene, and prints the reaction details.
    """
    # Define reactants using SMILES (Simplified Molecular Input Line Entry System) strings
    # Butadiene: C=CC=C
    # 1,1-dichloro-2,2-difluoroethene: ClC(Cl)=C(F)F
    reactant_smiles = ["C=CC=C", "ClC(Cl)=C(F)F"]
    reactant_names = ["Butadiene", "1,1-dichloro-2,2-difluoroethene"]

    # Define the generic Diels-Alder reaction using a SMARTS string
    # It maps how atoms in the diene and dienophile connect to form the cyclohexene product
    reaction_smarts = "[*:1]=[*:2][*:3]=[*:4].[*:5]=[*:6]>>[*:1]1[*:2]=[*:3][*:4][*:5][*:6]1"
    
    # Create RDKit molecule objects from SMILES
    reactants = [Chem.MolFromSmiles(s) for s in reactant_smiles]

    # Create and run the reaction
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    products = rxn.RunReactants(reactants)

    if products:
        # We expect one set of products, and one product in that set
        product_mol = products[0][0]
        Chem.SanitizeMol(product_mol) # Standardize the molecular representation
        
        # Calculate molecular formulas
        reactant_formulas = [CalcMolFormula(mol) for mol in reactants]
        product_formula = CalcMolFormula(product_mol)
        product_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"

        print(f"The product of the reaction between {reactant_names[0]} and {reactant_names[1]} is {product_name}.\n")

        # Construct and print the final equation
        final_equation = f"{reactant_formulas[0]} + {reactant_formulas[1]} -> {product_formula}"
        print("The chemical equation based on molecular formulas is:")
        print(final_equation)
        print()

        # Extract and print each number from the equation string
        print("The numbers found in the final equation are:")
        numbers = re.findall(r'\d', final_equation)
        for num in numbers:
            print(num)
    else:
        print("Could not determine the reaction product.")

if __name__ == "__main__":
    get_diels_alder_product()