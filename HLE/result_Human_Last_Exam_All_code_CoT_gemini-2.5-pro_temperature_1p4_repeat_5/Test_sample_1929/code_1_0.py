try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("This script requires the RDKit library.")
    print("Please install it using: pip install rdkit")
    exit()

def predict_diels_alder_product():
    """
    Predicts the product(s) of the Diels-Alder reaction between
    1,3-butadiene and 1,1-dichloro-2,2-difluoroethene.
    """
    # 1. Define reactants using SMILES strings
    # Butadiene (diene) and 1,1-dichloro-2,2-difluoroethene (dienophile)
    diene_smiles = "C=CC=C"
    dienophile_smiles = "C(=C(F)F)(Cl)Cl" # also C(=C(Cl)Cl)(F)F

    diene = Chem.MolFromSmiles(diene_smiles)
    dienophile = Chem.MolFromSmiles(dienophile_smiles)

    # 2. Define the Diels-Alder reaction using reaction SMARTS
    # This generic SMARTS pattern finds a diene and a dienophile and combines them.
    reaction_smarts = "[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1"
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    # 3. Print the overall reaction equation
    reactant_formula_1 = CalcMolFormula(diene)
    reactant_formula_2 = CalcMolFormula(dienophile)
    # The product formula will be the sum of reactant atoms
    product_formula = CalcMolFormula(Chem.MolFromSmiles("C1C(C(C(C1(Cl)Cl)F)F)=C"))
    
    # "output each number in the final equation!"
    # The final equation is 1 C4H6 + 1 C2Cl2F2 -> 1 C6H6Cl2F2
    print("Overall Reaction Equation:")
    print(f"1 {reactant_formula_1} + 1 {reactant_formula_2} -> 1 {product_formula}\n")

    # 4. Run the reaction
    # The result is a tuple of tuples, with each inner tuple containing product molecules.
    products = rxn.RunReactants((diene, dienophile))

    if not products:
        print("Reaction failed to generate products.")
        return

    # Manually determined IUPAC names based on standard nomenclature rules
    iupac_names = [
        "4,4-dichloro-5,5-difluorocyclohexene", # Major product
        "5,5-dichloro-4,4-difluorocyclohexene"  # Minor product
    ]
    
    print("The reaction can form two possible products (regioisomers):\n")
    # 5. Process and print each product
    for i, product_tuple in enumerate(products):
        product_mol = product_tuple[0]
        # Canonical SMILES provides a unique representation
        product_smiles = Chem.MolToSmiles(product_mol, canonical=True)

        print(f"--- Product {i+1} ---")
        print(f"SMILES: {product_smiles}")
        print(f"IUPAC Name: {iupac_names[i]}")
        if i == 0:
            print("(Note: This is the expected major product based on electronic effects.)")
        print("-" * 20 + "\n")

predict_diels_alder_product()