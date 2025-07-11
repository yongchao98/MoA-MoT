import sys

def generate_smiles():
    """
    This function generates and prints the SMILES representation of the molecule.
    The molecule is constructed based on the primary constraints:
    - Molecular Formula: C8H18N6
    - Core Functional Groups: 1 Azo group, 2 Amidine groups
    - Acyclic structure
    
    The resulting structure is 1,1'-(E)-diazene-1,2-diylbis(N',N'-dimethylmethanimidamide).
    SMILES string represents this structure.
    """
    # SMILES representation of the molecule
    smiles = "CN(C)C(=N)CN=NCC(=N)N(C)C"
    
    # Print the final SMILES string
    print(smiles)

# The following block allows the script to be run from the command line.
if __name__ == "__main__":
    generate_smiles()
