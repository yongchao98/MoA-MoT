import sys

def solve_chemistry_problem():
    """
    Solves a multi-step organic synthesis problem using RDKit.
    The reaction sequence is:
    1. Terpinolene + m-CPBA -> Compound 1 (Epoxide)
    2. Compound 1 + N,N-dimethyl thioformamide -> Compound 2 (Thiirane)
    3. Compound 2 + LiAlH4 -> Compound 3 (Thiol)
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import Descriptors
    except ImportError:
        print("Error: RDKit library not found.")
        print("Please install it to run this script, for example, using pip:")
        print("pip install rdkit")
        sys.exit(1)

    # Step 0: Define the starting material, terpinolene
    # IUPAC Name: 4-isopropylidene-1-methylcyclohex-1-ene
    terpinolene_smiles = 'CC1=CCC(CC1)=C(C)C'
    terpinolene = Chem.MolFromSmiles(terpinolene_smiles)
    terpinolene_formula = Descriptors.CalcMolFormula(terpinolene)
    
    print("--- Reaction Sequence Analysis ---")
    print(f"\nStep 0: Starting Material")
    print(f"Name: Terpinolene")
    print(f"Formula: {terpinolene_formula}")
    print(f"SMILES: {terpinolene_smiles}")
    print("-" * 35)

    # Step 1: Reaction of terpinolene with m-CPBA (Epoxidation)
    # m-CPBA selectively epoxidizes the endocyclic (in-ring) double bond of terpinolene.
    rxn1_smarts = '[C;R:1]=[C;R:2]>>[C:1]1O[C:2]1'
    rxn1 = AllChem.ReactionFromSmarts(rxn1_smarts)
    products1 = rxn1.RunReactants((terpinolene,))
    compound1 = products1[0][0]
    Chem.SanitizeMol(compound1)
    compound1_smiles = Chem.MolToSmiles(compound1)
    compound1_formula = Descriptors.CalcMolFormula(compound1)

    print(f"Step 1: Epoxidation with m-CPBA")
    print(f"Product (Compound 1): Terpinolene epoxide")
    print(f"Formula: {compound1_formula}")
    print(f"SMILES: {compound1_smiles}")
    print("-" * 35)

    # Step 2: Reaction of Compound 1 to form Compound 2 (Thiirane)
    # The epoxide is converted to a thiirane (episulfide).
    rxn2_smarts = '[C:1]1O[C:2]1>>[C:1]1S[C:2]1'
    rxn2 = AllChem.ReactionFromSmarts(rxn2_smarts)
    products2 = rxn2.RunReactants((compound1,))
    compound2 = products2[0][0]
    Chem.SanitizeMol(compound2)
    compound2_smiles = Chem.MolToSmiles(compound2)
    compound2_formula = Descriptors.CalcMolFormula(compound2)

    print(f"Step 2: Conversion to Thiirane")
    print(f"Product (Compound 2): Terpinolene thiirane")
    print(f"Formula: {compound2_formula}")
    print(f"SMILES: {compound2_smiles}")
    print("-" * 35)

    # Step 3: Reduction of Compound 2 with LiAlH4
    # LiAlH4 reduces the thiirane to a thiol. The hydride attacks the less substituted
    # carbon of the thiirane ring, leading to the formation of a thiol at the
    # more substituted carbon.
    # SMARTS for this regioselective ring-opening:
    # [C;X4;H0:1] = tertiary carbon of the thiirane
    # [C;X3;H1:2] = secondary carbon of the thiirane
    rxn3_smarts = '[C;X4;H0:1]1S[C;X3;H1:2]1>>[C:1](S)[C:2]'
    rxn3 = AllChem.ReactionFromSmarts(rxn3_smarts)
    products3 = rxn3.RunReactants((compound2,))
    
    if not products3:
        print("Error: Could not determine the product of Step 3.")
        sys.exit(1)

    compound3 = products3[0][0]
    Chem.SanitizeMol(compound3)
    compound3_smiles = Chem.MolToSmiles(compound3)
    compound3_formula = Descriptors.CalcMolFormula(compound3)
    compound3_mw = Descriptors.ExactMolWt(compound3)
    
    print(f"Step 3: Reduction with LiAlH4")
    print(f"Final Product (Compound 3)")
    print(f"Name: 1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol")
    print(f"Formula: {compound3_formula}")
    print(f"Molecular Weight: {compound3_mw:.4f}")
    print(f"SMILES: {compound3_smiles}")
    print("-" * 35)
    
    # The request was to "output each number in the final equation"
    # This can be interpreted as showing the change in molecular formula.
    print("\n--- Final Answer Summary ---")
    print("The final compound, Compound 3, is 1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol.")
    print("\nThe reaction equation in terms of molecular formulas is:")
    final_equation = f"{terpinolene_formula}  ->  {compound1_formula}  ->  {compound2_formula}  ->  {compound3_formula}"
    print(final_equation)


if __name__ == '__main__':
    solve_chemistry_problem()
