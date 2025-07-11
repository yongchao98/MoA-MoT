import re
# The rdkit library is a standard toolkit for cheminformatics.
# If you don't have it, you can install it via pip: pip install rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def solve_wittig_reaction():
    """
    This function performs a Wittig reaction between pivalaldehyde and a
    specified phosphonium ylide, and prints the results.
    """
    # Step 1: Define the reactants using their SMILES string representations.
    # Pivalaldehyde: (CH3)3C-CHO
    pivalaldehyde_smiles = "CC(C)(C)C=O"
    pivalaldehyde_name = "pivalaldehyde"

    # Ylide: (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane
    # Structure: Ph3P=CH-CH2-(ortho-chlorophenyl)
    ylide_smiles = "c1(Cl)ccccc1CC=[P](c2ccccc2)(c3ccccc3)c4ccccc4"
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    
    reactants = (Chem.MolFromSmiles(pivalaldehyde_smiles), Chem.MolFromSmiles(ylide_smiles))

    # Step 2: Define the Wittig reaction mechanism using a generic reaction pattern (SMARTS).
    # This pattern says: find a C=O and a C=P group, and swap the O and C groups.
    # [C:1]=[O:2] is the aldehyde. [#6:4]=[P:3] is the ylide carbon and phosphorus.
    # The ">>" indicates the transformation to the products.
    rxn_smarts = "[C:1]=[O:2].[#6:4]=[P:3]>>[C:1]=[#6:4].[O:2]=[P:3]"
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)

    # Step 3: Run the reaction simulation.
    product_sets = rxn.RunReactants(reactants)

    # The reaction yields the alkene and triphenylphosphine oxide. We need to identify them.
    alkene_product = None
    if product_sets:
        # We find the main organic product by filtering out the one with phosphorus.
        for mol in product_sets[0]:
            if not mol.HasSubstructMatch(Chem.MolFromSmarts("[P]")):
                alkene_product = mol
                break
    
    # Step 4: Determine the IUPAC name and present the results.
    # The structure formed is (CH3)3C-CH=CH-CH2-(2-Cl-Ph).
    # The systematic IUPAC name is 1-(2-chlorophenyl)-4,4-dimethylpent-2-ene.
    # Note: Stereochemistry (E/Z) is not specified as it depends on conditions.
    product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct_name = "triphenylphosphine oxide"

    print("--- Wittig Reaction Analysis ---")
    print(f"Reactant 1: {pivalaldehyde_name}")
    print(f"Reactant 2: {ylide_name}\n")
    print("The reaction equation is:")
    print(f"{pivalaldehyde_name} + {ylide_name} -> {product_name} + {byproduct_name}\n")
    
    print(f"The main organic product is named: {product_name}\n")
    
    # Final step: As requested, output each number from the product's IUPAC name.
    print("The numbers in the final product's IUPAC name equation are:")
    numbers = re.findall(r'\d+', product_name)
    for num in numbers:
        print(num)

solve_wittig_reaction()