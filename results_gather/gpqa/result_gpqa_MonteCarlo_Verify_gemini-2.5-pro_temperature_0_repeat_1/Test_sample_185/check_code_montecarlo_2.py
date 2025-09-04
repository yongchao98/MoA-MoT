from rdkit import Chem
from rdkit.Chem import AllChem

def check_cope_rearrangement():
    """
    Checks if the Cope rearrangement of (1S,4R)-2-vinyl-2-azabicyclo[2.2.1]hept-5-ene
    yields the product described in option A.
    """
    # 1. Define the reactant molecule using its SMILES string.
    # The SMILES for 2-vinyl-2-azabicyclo[2.2.1]hept-5-ene is C=CN1C2CC(C=C2)C1.
    # Stereochemistry (1S,4R) is ignored as we are checking for constitutional isomerism (connectivity).
    reactant_smiles = "C=CN1C2CC(C=C2)C1"
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    if not reactant_mol:
        return "Incorrect: Failed to parse reactant SMILES."

    # 2. Define the aza-Cope [3,3]-sigmatropic rearrangement using reaction SMARTS.
    # The pattern identifies the C=C-N-C-C=C system and defines the bond changes.
    # [C:1]=[C:2]-[N:3]-[C:4]-[C:5]=[C:6] >> [C:2]=[N:3]-[C:4]=[C:5]-[C:6]-[C:1]
    # Note: RDKit's SMARTS requires atom mapping to define the transformation.
    # The bond between atoms 3 and 4 is broken. A new bond between 1 and 6 is formed.
    # The double bonds shift.
    reaction_smarts = "[C:1]=[C:2][N:3][C:4][C:5]=[C:6]>>[C:1][C:6]=[C:5][C:4][N:3]=[C:2]"
    reaction = AllChem.ReactionFromSmarts(reaction_smarts)

    # 3. Run the reaction on the reactant.
    products = reaction.RunReactants((reactant_mol,))

    if not products:
        return "Incorrect: The defined Cope rearrangement did not apply to the reactant."

    # The reaction should yield one major product topology.
    if len(products) > 1:
        # This might happen if the molecule has multiple reacting sites, but not in this case.
        # We'll check the first product.
        pass

    product_mol = products[0][0]

    # 4. Sanitize the product molecule to correct valencies and aromatization flags.
    try:
        Chem.SanitizeMol(product_mol)
    except Exception as e:
        return f"Incorrect: Produced an invalid chemical structure. Sanitization error: {e}"

    # Generate the canonical SMILES for the simulated product.
    # Canonical SMILES provides a unique string for a given molecular graph,
    # making it perfect for isomorphism checking.
    product_canonical_smiles = Chem.MolToSmiles(product_mol, canonical=True)

    # 5. Define the expected product (Option A) and compare.
    # The product of the Cope rearrangement is a hexahydrocyclopenta[c]pyridine derivative.
    # Based on chemical principles, its structure corresponds to the SMILES 'C1=CC2C3CN=CC3C12'.
    # We assert that this is the structure represented by the name in Option A.
    option_A_smiles = "C1=CC2C3CN=CC3C12"
    option_A_mol = Chem.MolFromSmiles(option_A_smiles)
    option_A_canonical_smiles = Chem.MolToSmiles(option_A_mol, canonical=True)

    # 6. Check for correctness.
    if product_canonical_smiles == option_A_canonical_smiles:
        return "Correct"
    else:
        return (f"Incorrect: The simulated product does not match the structure of Option A.\n"
                f"Simulated Product SMILES: {product_canonical_smiles}\n"
                f"Expected Option A SMILES: {option_A_canonical_smiles}")

# Run the check
result = check_cope_rearrangement()
print(result)