# To run this code, you need to install the RDKit library:
# pip install rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def get_formula(smiles):
    """Calculates the molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES"
    return CalcMolFormula(mol)

def main():
    """
    Solves the chemical reaction problem by identifying Compound A and
    displaying the balanced chemical equation.
    """
    # SMILES strings for reactants and products
    smiles_aminopyridine = "Nc1ncccc1"
    smiles_opa = "O=Cc1ccccc1C=O"
    smiles_hcn = "C#N"
    smiles_compound_A = "n1(c2ccccc2c1C#N)c3ncccc3" # 2-(pyridin-2-yl)isoindole-1-carbonitrile
    smiles_water = "O"

    # Get molecular formulas
    formula_aminopyridine = get_formula(smiles_aminopyridine)
    formula_opa = get_formula(smiles_opa)
    formula_hcn = get_formula(smiles_hcn)
    formula_compound_A = get_formula(smiles_compound_A)
    formula_water = get_formula(smiles_water)

    # --- Outputting the results ---
    print("The reaction is a multi-component reaction leading to a heterocyclic product.")
    print("The overall balanced chemical equation is:\n")
    
    # Printing the equation with coefficients
    print(f"1 {formula_aminopyridine} (2-Aminopyridine) + 1 {formula_opa} (o-Phthalaldehyde) + 1 {formula_hcn} (Cyanide Source) -> 1 {formula_compound_A} (Compound A) + 2 {formula_water} (Water)")

    print("\n----------------------------------------")
    print("Conclusion: What is Compound A?")
    print("----------------------------------------")
    print(f"Compound A is 2-(pyridin-2-yl)isoindole-1-carbonitrile.")
    print(f"Its molecular formula is: {formula_compound_A}")
    print(f"Its SMILES string is: {smiles_compound_A}")

if __name__ == "__main__":
    main()
