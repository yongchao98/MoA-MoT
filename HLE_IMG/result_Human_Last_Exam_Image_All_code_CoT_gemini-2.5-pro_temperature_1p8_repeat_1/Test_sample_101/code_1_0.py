# To run this code, you may need to install the RDKit library.
# You can do this by running: pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, GetAtomCounts

def identify_compound_A():
    """
    Identifies the final product of a two-step organic synthesis and calculates its properties.
    """

    print("Step-by-step analysis of the reaction:")
    print("1. The first step is the condensation of 3-hydroxy-pyridine-2-carbaldehyde and aniline. This forms an imine intermediate.")
    print("2. The second step is a Strecker-type reaction where sodium cyanide (NaCN) is added to the imine.")
    print("   The cyanide ion (CN-) attacks the carbon of the imine's C=N bond, and the nitrogen is subsequently protonated.")
    print("-" * 50)
    
    # The structure of compound A is (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile.
    # We can represent this molecule using a SMILES string.
    product_smiles = "N#CC(Nc1ccccc1)c1ncccc1O"
    product_name = "(3-hydroxypyridin-2-yl)(phenylamino)acetonitrile"
    
    print(f"The final product, Compound A, is identified as: {product_name}")
    print(f"SMILES string: {product_smiles}\n")
    
    # Create an RDKit molecule object from the SMILES string
    mol = Chem.MolFromSmiles(product_smiles)
    
    if mol:
        # Calculate molecular properties
        molecular_formula = CalcMolFormula(mol)
        exact_mass = Descriptors.ExactMolWt(mol)
        
        print("Chemical Properties of Compound A:")
        print(f"Molecular Formula: {molecular_formula}")
        print("\nThe numbers in the final chemical formula are:")
        
        # GetAtomCounts returns a dictionary of atomic number to count.
        # We convert atomic numbers to symbols for clarity.
        atom_counts = GetAtomCounts(mol)
        for atomic_num in sorted(atom_counts.keys()):
            atom_symbol = Chem.Atom(atomic_num).GetSymbol()
            count = atom_counts[atomic_num]
            print(f"  - Number of {atom_symbol} atoms: {count}")
        
        print(f"\nExact Molar Mass: {exact_mass:.4f} g/mol")
        
    else:
        print("Error: Could not generate the molecule from the SMILES string.")

# Run the function to get the answer
identify_compound_A()