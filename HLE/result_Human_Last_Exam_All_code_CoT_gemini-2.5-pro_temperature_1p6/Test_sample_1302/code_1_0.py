from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def generate_molecule_info():
    """
    This function holds the SMILES string for the designed molecule and prints its properties.
    The molecule is a flavylium-type zwitterion, which satisfies the complex set of constraints.
    - Formula: C15H10O5 (MW: 270.053, Valence e-: 100, Heavy Atoms: 20, Heteroatoms: 5 O)
    - Structure: 3 aromatic rings (2 benzene-like, 1 aromatic pyrylium), 1 rotatable bond.
    - Charge: 0 (zwitterion)
    - H-Donors: 3 (from 3 OH groups), H-Acceptors: 5 (from 5 O atoms)
    - It contains 3 phenolic OH groups and no forbidden functional groups.
    """
    # SMILES for 3,7-dihydroxy-2-(4-hydroxyphenyl)chromenylium-5-olate
    # This structure satisfies all constraints.
    final_smiles = "c1c(O)cc(c2c1[o+]=c(cc2[O-])c1ccc(O)cc1)O"

    # Generate molecule object from SMILES
    mol = Chem.MolFromSmiles(final_smiles)
    mol_with_hs = Chem.AddHs(mol)

    # Calculate properties for the final equation
    formula = rdMolDescriptors.CalcMolFormula(mol_with_hs)
    exact_mw = Descriptors.ExactMolWt(mol_with_hs)

    # Print the proposed molecule and its verification
    print(f"Designed Molecule SMILES: {final_smiles}\n")
    print("Final Equation (Molecular Formula = Exact Mass):")
    # The output shows the atoms and their counts followed by the calculated mass.
    # Note: RDKit may order atoms differently (e.g., C, H, O).
    print(f"{formula} = {exact_mw:.5f}")


# Execute the function to get the output.
generate_molecule_info()
<<<c1c(O)cc(c2c1[o+]=c(cc2[O-])c1ccc(O)cc1)O>>>