try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
except ImportError:
    print("Error: RDKit library not found.")
    print("Please install it to run this code, for example using pip:")
    print("pip install rdkit")
    exit()

def solve_reaction():
    """
    Solves the given two-step chemical reaction to identify Compound A.
    """
    print("Analyzing the chemical reaction to determine Compound A...")
    print("-" * 50)
    
    # Step 1: Imine Formation
    print("Step 1: Reaction of 3-hydroxy-pyridine-2-carbaldehyde with aniline.")
    print("This is a condensation reaction forming an imine intermediate.")
    reactant1_smiles = "O=Cc1c(O)cccn1"  # 3-hydroxy-pyridine-2-carbaldehyde
    reactant2_smiles = "Nc1ccccc1"      # Aniline
    intermediate_smiles = "Oc1ncccc1C=Nc1ccccc1" # (E)-2-((phenylimino)methyl)pyridin-3-ol
    
    print(f"Reactant 1 (SMILES): {reactant1_smiles}")
    print(f"Reactant 2 (SMILES): {reactant2_smiles}")
    print(f"Intermediate (SMILES): {intermediate_smiles}")
    print("-" * 50)

    # Step 2: Cyanide Addition to the Imine
    print("Step 2: Reaction of the imine intermediate with NaCN.")
    print("This is a nucleophilic addition of cyanide (CN-) to the C=N bond.")
    print("The result is an alpha-aminonitrile.")
    
    # Create an RDKit molecule object for the intermediate
    intermediate_mol = Chem.MolFromSmiles(intermediate_smiles)
    if intermediate_mol is None:
        print("Error: Invalid intermediate SMILES string.")
        return

    # Define the reaction for cyanide addition using SMARTS.
    # This pattern finds a C=N bond and adds a cyano group to the C
    # and a hydrogen to the N.
    # [C:1]=[N:2] >> [C:1]([C]#[N])[N:2H1]
    rxn_smarts = '[C:1]=[N:2]>>[C:1]([C]#[N])[N:2H1]'
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)

    # Run the reaction
    products = rxn.RunReactants((intermediate_mol,))

    if not products:
        print("Error: The reaction simulation failed.")
        return
        
    # Process the product
    # The result is a tuple of tuples; we want the first molecule of the first result set.
    product_mol = products[0][0]
    Chem.SanitizeMol(product_mol)

    # Get the canonical SMILES and molecular formula for the final product
    final_smiles = Chem.MolToSmiles(product_mol, canonical=True)
    final_formula = rdMolDescriptors.CalcMolFormula(product_mol)

    print("-" * 50)
    print("Result: Compound A has been identified.")
    print("\n--- Compound A ---")
    print(f"Chemical Structure (SMILES): {final_smiles}")
    print(f"Molecular Formula: {final_formula}")
    print("IUPAC Name: 2-((3-hydroxypyridin-2-yl)(phenylamino))acetonitrile")
    print("\nBreakdown of the molecular formula calculation:")
    # Deconstruct the formula string to print each number as requested
    import re
    atom_counts = re.findall(r'([A-Z][a-z]*)(\d*)', final_formula)
    equation_parts = []
    for atom, count in atom_counts:
        num = count if count else '1'
        equation_parts.append(f"{atom}: {num}")
    print(" + ".join(equation_parts))


if __name__ == '__main__':
    solve_reaction()
