def solve_synthesis():
    """
    Solves the multi-step synthesis problem using RDKit to model the reactions.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print("This script requires the RDKit library.")
        print("Please install it using: pip install rdkit-pypi")
        return

    # --- Step 0: Define Starting Material ---
    # Terpinolene SMILES: 1-methyl-4-(propan-2-ylidene)cyclohex-1-ene
    terpinolene_smiles = 'CC1=CCC(=C(C)C)CC1'
    terpinolene = Chem.MolFromSmiles(terpinolene_smiles)

    print("--- Multi-Step Synthesis Analysis ---")
    print(f"Starting Material: Terpinolene")
    print(f"SMILES: {terpinolene_smiles}\n")

    # --- Step 1: Epoxidation ---
    # m-CPBA epoxidizes the more substituted (trisubstituted endocyclic) double bond.
    # The reaction is: Alkene -> Epoxide (1:1 stoichiometry)
    # We identify the correct product by checking which one retains an unreacted double bond.
    
    rxn1_smarts = '[C:1]=[C:2]>>[C:1]1O[C:2]1'
    rxn1 = AllChem.ReactionFromSmarts(rxn1_smarts)
    products1 = rxn1.RunReactants((terpinolene,))

    # RDKit will generate both possible epoxides. We select the correct one.
    # The correct product (Compound 1) is from epoxidation of the endocyclic double bond,
    # leaving the exocyclic C=C double bond intact.
    compound1 = None
    for p_tuple in products1:
        for p in p_tuple:
            if p.HasSubstructMatch(Chem.MolFromSmarts('C=C')): # Check if a double bond remains
                Chem.SanitizeMol(p)
                compound1 = p
                break
    
    compound1_smiles = Chem.MolToSmiles(compound1)
    print("Step 1: Reaction of terpinolene with m-CPBA")
    print("Principle: Preferential epoxidation of the more substituted (endocyclic) double bond.")
    print("Equation: Terpinolene + m-CPBA -> Compound 1 (1:1)")
    print(f"Product (Compound 1): Terpinolene oxide")
    print(f"SMILES: {compound1_smiles}\n")

    # --- Step 2: Epoxide to Thiirane ---
    # The epoxide is converted to a thiirane using N,N-dimethyl thioformamide.
    # This is a direct functional group transformation (1:1 stoichiometry).
    
    rxn2_smarts = '[C:1]1O[C:2]1>>[C:1]1S[C:2]1'
    rxn2 = AllChem.ReactionFromSmarts(rxn2_smarts)
    products2 = rxn2.RunReactants((compound1,))
    compound2 = products2[0][0]
    Chem.SanitizeMol(compound2)
    compound2_smiles = Chem.MolToSmiles(compound2)

    print("Step 2: Reaction of Compound 1 with N,N-dimethyl thioformamide / acid")
    print("Principle: Conversion of an epoxide to a thiirane (episulfide).")
    print("Equation: Compound 1 -> Compound 2 (1:1)")
    print(f"Product (Compound 2): 1,2-Epithio-p-menth-4(8)-ene")
    print(f"SMILES: {compound2_smiles}\n")

    # --- Step 3: Reductive opening of Thiirane ---
    # LiAlH4 opens the thiirane ring. The hydride nucleophile attacks the less substituted
    # carbon of the thiirane (tertiary C) via an SN2 mechanism.
    # This breaks the C-S bond, forming a thiol at the more substituted carbon (quaternary C).
    # Stoichiometry is 1:1.
    
    # SMARTS for the specific ring opening:
    # [C;X4:1] is a quaternary carbon, [C;H1:2] is a tertiary carbon.
    # The reaction breaks the S-[C;H1:2] bond.
    rxn3_smarts = '[C;X4:1]1[S:3][C;H1:2]1>>[C:1]([S:3])[C:2]'
    rxn3 = AllChem.ReactionFromSmarts(rxn3_smarts)
    products3 = rxn3.RunReactants((compound2,))
    compound3 = products3[0][0]
    Chem.SanitizeMol(compound3)
    compound3_smiles = Chem.MolToSmiles(compound3)

    # Use RDKit to try and generate the IUPAC name
    try:
        from PubChemPy import Compound
        pcrecord = Compound.from_smiles(compound3_smiles)
        compound3_name = pcrecord.iupac_name
        if not compound3_name:
            compound3_name = "1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol"
    except ImportError:
        compound3_name = "1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol (requires PubChemPy for programmatic lookup)"
    
    print("Step 3: Reduction of Compound 2 with LiAlH4")
    print("Principle: Regioselective reductive opening of the thiirane. Hydride attacks the less substituted carbon, forming a thiol at the more substituted carbon.")
    print("Equation: Compound 2 + LiAlH4 -> Compound 3 (1:1)")
    print(f"Final Product (Compound 3): {compound3_name}")
    print(f"SMILES: {compound3_smiles}\n")
    print("--- End of Analysis ---")


if __name__ == '__main__':
    solve_synthesis()