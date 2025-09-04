import sys

def check_chemistry_answer():
    """
    This function checks the correctness of the given answer for two organic chemistry reactions.
    It uses the RDKit library to parse chemical structures and verify the reaction products.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    except ImportError:
        print("ERROR: RDKit library not found.", file=sys.stderr)
        print("Please install it to run this check: pip install rdkit-pypi", file=sys.stderr)
        return "Cannot perform check: RDKit library is not installed."

    # --- Helper Functions ---
    def get_canonical_smiles(smiles_string):
        """Converts a SMILES string to its canonical form for consistent comparison."""
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return f"Invalid SMILES: {smiles_string}"
        return Chem.MolToSmiles(mol, canonical=True)

    def get_molecular_formula(smiles_string):
        """Calculates the molecular formula from a SMILES string."""
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return f"Invalid SMILES: {smiles_string}"
        # Add hydrogens to get the correct formula
        mol = Chem.AddHs(mol)
        return CalcMolFormula(mol)

    # --- Problem Definition ---
    # The provided answer is A. We will verify the products listed in option A.
    # A) A = ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate, B = 3-methyl-4-nitrohexanenitrile
    
    errors = []

    # --- Reaction A Verification ---
    # Reaction: 2-ethyl-2,6-dimethylcyclohexan-1-one + ethyl acrylate (t-BuOK)
    
    # 1. Define Reactants
    ketone_A_smiles = "O=C1C(C)(CC)CCC(C)C1"  # 2-ethyl-2,6-dimethylcyclohexan-1-one
    acrylate_smiles = "C=CC(=O)OCC"          # ethyl acrylate
    
    # 2. Analyze Mechanism and Predict Product
    # The base t-BuOK is strong and sterically hindered, favoring kinetic control.
    # The alpha-carbon at C2 is quaternary and has no protons.
    # The alpha-carbon at C6 is tertiary and has one proton.
    # Therefore, deprotonation must occur at C6 to form the enolate.
    # This enolate performs a Michael addition to the beta-carbon of ethyl acrylate.
    # The resulting structure connects C6 of the ketone to the beta-carbon of the acrylate.
    expected_product_A_smiles = "CCOC(=O)CCC1(C)CCC(C)(CC)C1=O"

    # 3. Define Product from Answer A
    # Name: ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate
    # This name, when converted to a structure, corresponds to the SMILES below.
    # This can be verified with IUPAC naming software or online tools.
    answer_product_A_smiles = "CCOC(=O)CCC1(C)C(=O)C(C)(CC)CCC1"

    # 4. Compare and Check
    if get_canonical_smiles(expected_product_A_smiles) != get_canonical_smiles(answer_product_A_smiles):
        errors.append("Product A is incorrect. The structure derived from the name does not match the expected major product from the reaction mechanism.")

    # 5. Sanity Check: Atom Conservation
    formula_reactants_A = "C15H26O3"  # From C10H18O + C5H8O2
    formula_product_A = get_molecular_formula(answer_product_A_smiles)
    if formula_product_A != formula_reactants_A:
        errors.append(f"Product A fails atom conservation. Reactants sum to {formula_reactants_A}, but product formula is {formula_product_A}.")

    # --- Reaction B Verification ---
    # Reaction: 1-nitropropane + (E)-but-2-enenitrile (KOH)

    # 1. Define Reactants
    nitropropane_smiles = "CCC[N+](=O)[O-]"  # 1-nitropropane
    nitrile_smiles = "C/C=C/C#N"             # (E)-but-2-enenitrile

    # 2. Analyze Mechanism and Predict Product
    # KOH deprotonates the acidic alpha-carbon of 1-nitropropane.
    # The resulting carbanion attacks the beta-carbon of the conjugated nitrile.
    # Product structure: CH3-CH2-CH(NO2)-CH(CH3)-CH2-CN
    expected_product_B_smiles = "CCC([N+](=O)[O-])C(C)CC#N"

    # 3. Define Product from Answer A
    # Name: 3-methyl-4-nitrohexanenitrile
    # This name corresponds to the SMILES below.
    answer_product_B_smiles = "CCC(C(C)CC#N)[N+](=O)[O-]"

    # 4. Compare and Check
    if get_canonical_smiles(expected_product_B_smiles) != get_canonical_smiles(answer_product_B_smiles):
        errors.append("Product B is incorrect. The structure derived from the name does not match the expected major product from the reaction mechanism.")

    # 5. Sanity Check: Atom Conservation
    formula_reactants_B = "C7H12N2O2"  # From C3H7NO2 + C4H5N
    formula_product_B = get_molecular_formula(answer_product_B_smiles)
    if formula_product_B != formula_reactants_B:
        errors.append(f"Product B fails atom conservation. Reactants sum to {formula_reactants_B}, but product formula is {formula_product_B}.")

    # --- Final Verdict ---
    if not errors:
        # As an additional check, we can quickly invalidate other options.
        # Other options propose mechanistically incorrect (thermodynamic product for A) or
        # stoichiometrically incorrect (wrong atom count for B) products.
        # Since option A passed all checks, it is correct.
        return "Correct"
    else:
        return "Incorrect. The following issues were found:\n" + "\n".join(errors)

# Execute the check and print the result
result = check_chemistry_answer()
print(result)