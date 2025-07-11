import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def get_molecule_info(mol):
    """Helper function to get molecule properties."""
    # Sanitize the molecule to ensure correct properties
    Chem.SanitizeMol(mol)
    smiles = Chem.MolToSmiles(mol)
    formula = CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    return smiles, formula, mw

def print_step(step_number, description, reactant_name, product_name, product_mol):
    """Helper function to print information for each reaction step."""
    smiles, formula, mw = get_molecule_info(product_mol)
    print(f"### Step {step_number}: {description} ###")
    print(f"Reactant: {reactant_name}")
    print(f"Product: {product_name}")
    print(f"  - SMILES: {smiles}")
    print(f"  - Formula: {formula}")
    print(f"  - Molecular Weight: {mw:.2f}\n")

def main():
    """
    Main function to run the reaction sequence analysis.
    """
    print("--- Analyzing the 3-Step Chemical Synthesis ---\n")

    # Step 0: Starting Material - Terpinolene
    terpinolene_smiles = "CC1=CCC(CC1)=C(C)C"
    terpinolene = Chem.MolFromSmiles(terpinolene_smiles)
    smiles, formula, mw = get_molecule_info(terpinolene)
    print("### Step 0: Starting Material ###")
    print("Reactant: Terpinolene")
    print(f"  - SMILES: {smiles}")
    print(f"  - Formula: {formula}")
    print(f"  - Molecular Weight: {mw:.2f}\n")

    # Step 1: Epoxidation with m-CPBA
    # m-CPBA selectively epoxidizes the more reactive endocyclic double bond.
    # The product is terpinolene oxide.
    compound1_smiles = "CC12OC1CCC(CC2)=C(C)C"
    compound1 = Chem.MolFromSmiles(compound1_smiles)
    print_step(1, "Epoxidation with m-CPBA", "Terpinolene", "Compound 1", compound1)

    # Step 2: Conversion of Epoxide to Thiirane
    # N,N-dimethyl thioformamide replaces the epoxide oxygen with sulfur.
    rxn_smarts_step2 = "[C:1]1O[C:2]1>>[C:1]1S[C:2]1"
    rxn2 = AllChem.ReactionFromSmarts(rxn_smarts_step2)
    products2 = rxn2.RunReactants((compound1,))
    compound2 = products2[0][0]
    print_step(2, "Thiirane formation", "Compound 1", "Compound 2", compound2)

    # Step 3: Reduction of Thiirane with LiAlH4
    # LiAlH4 opens the thiirane ring to form a thiol. The hydride attacks the
    # less hindered carbon, forming the thiol at the more substituted carbon.
    # The product is 1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol.
    compound3_smiles = "CC(S)(C1)CCC(=C(C)C)CC1"
    compound3 = Chem.MolFromSmiles(compound3_smiles)
    print_step(3, "Reduction with LiAlH4", "Compound 2", "Compound 3", compound3)

    # --- Final Answer ---
    print("--- Conclusion ---")
    final_name = "1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol"
    print(f"The final product, Compound 3, is: {final_name}")
    print("This is also known as p-menth-4(8)-ene-1-thiol or 'grapefruit mercaptan'.\n")

    # As requested, outputting the numbers from the final product's chemical formula.
    final_formula = CalcMolFormula(compound3)
    print(f"The molecular formula of the final product is {final_formula}.")
    
    # Parse the formula string e.g., "C10H18S" into its components
    parts = re.findall(r'([A-Z][a-z]*)(\d*)', final_formula)
    equation_parts = []
    for element, count in parts:
        # If count is empty, it means 1 atom
        num = int(count) if count else 1
        equation_parts.append(f"{element} = {num}")
    
    print("The final equation for the molecular formula is:")
    print(" + ".join(equation_parts))


if __name__ == "__main__":
    main()