# The user might need to install rdkit first, e.g., via pip: pip install rdkit
import re
try:
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("This script requires the RDKit library.")
    print("Please install it, for example using: pip install rdkit")
else:
    # The reaction is a Joullié-Ugi three-component reaction.
    # The product, compound A, is 1-cyano-2-(pyridin-2-yl)isoindolin-3-ol.
    # We can represent this compound using its SMILES string.
    # The constitution is represented, ignoring stereoisomers which would be formed as a mixture.
    smiles_A = "N#CC1c2ccccc2C(O)N1c1ncccc1"

    # Create an RDKit molecule object from the SMILES string
    molecule_A = Chem.MolFromSmiles(smiles_A)

    if molecule_A:
        # Calculate the molecular formula
        molecular_formula = CalcMolFormula(molecule_A)

        print("The reaction between 2-aminopyridine, o-phthalaldehyde, and TMSCN produces compound A.")
        print("Compound A is: 1-cyano-2-(pyridin-2-yl)isoindolin-3-ol")
        print(f"Molecular Formula: {molecular_formula}")

        # The reaction conditions given are: 28 C, 4 h
        print("Reaction conditions: Temperature = 28 C, Time = 4 h.")
        
        # As requested, outputting the numbers from the final formula "equation"
        # For the formula C14H11N3O, the numbers are the atom counts: 14, 11, 3, and 1.
        print("The numbers in the molecular formula are:")
        
        # Use regex to find all element symbols and their counts
        atom_counts = re.findall(r'([A-Z][a-z]*)(\d*)', molecular_formula)
        
        final_numbers = []
        for element, count in atom_counts:
            # If count is empty, it means the count is 1
            if not count:
                final_numbers.append('1')
            else:
                final_numbers.append(count)
        
        print(' '.join(final_numbers))

    else:
        print("Error: Could not process the chemical structure.")