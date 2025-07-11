try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.Draw import MolToImage
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit is not installed. Will proceed with text-based output.")
    print("To visualize molecules, please install rdkit: pip install rdkit\n")

def get_product_info(name, smiles):
    """Generates a description for a chemical product."""
    info = f"Name: {name}\nSMILES: {smiles}"
    if RDKIT_AVAILABLE:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
            mw = Descriptors.ExactMolWt(mol)
            info += f"\nMolecular Formula: {formula}\nExact Mass: {mw:.4f}"
        else:
            info += "\nCould not generate molecule from SMILES."
    return info

def solve_reaction():
    """
    Identifies and describes the products of the given chemical reaction.
    """
    # The two major products are constitutional isomers.
    # We will label them A and B. The assignment is arbitrary as the problem
    # does not provide data to distinguish them (e.g., GC retention time).

    # Product A: 2-(tert-butoxy)-1-phenylethyl benzoate
    # Formed by adding -OtBu to CH2 and -OCOPh to the benzylic CH.
    smiles_A = "c1ccccc1C(OC(=O)c2ccccc2)COC(C)(C)C"
    name_A = "2-(tert-butoxy)-1-phenylethyl benzoate"

    # Product B: 1-(tert-butoxy)-1-phenylethyl benzoate
    # Formed by adding -OCOPh to CH2 and -OtBu to the benzylic CH.
    smiles_B = "c1ccccc1C(OC(C)(C)C)COC(=O)c2ccccc2"
    name_B = "1-(tert-butoxy)-1-phenylethyl benzoate"
    
    print("The reaction of styrene with tert-butyl peroxybenzoate results in the 1,2-addition of a tert-butoxy group and a benzoyloxy group across the double bond.")
    print("This yields two major products, A and B, which are constitutional isomers.\n")

    print("--- Product A ---")
    print(get_product_info(name_A, smiles_A))
    print("\n-------------------\n")

    print("--- Product B ---")
    print(get_product_info(name_B, smiles_B))
    print("\n-------------------\n")

if __name__ == "__main__":
    solve_reaction()