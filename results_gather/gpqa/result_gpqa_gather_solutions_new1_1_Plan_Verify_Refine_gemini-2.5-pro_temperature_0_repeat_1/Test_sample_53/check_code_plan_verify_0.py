import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def check_answer():
    """
    Checks the correctness of the answer for the RCM synthesis question.
    """
    try:
        # 1. Define the target product and its canonical SMILES representation.
        # The name is 5-isopropyl-3,4-dimethylcyclohex-1-ene
        # A valid SMILES for this is C1=CC(C)C(C)C(C(C)C)C1
        target_mol = Chem.MolFromSmiles("C1=CC(C)C(C)C(C(C)C)C1")
        if not target_mol:
            return "Error: Could not parse the target molecule's structure."
        target_canonical_smiles = Chem.MolToSmiles(target_mol, canonical=True)

        # 2. Define the candidate starting materials from the options.
        options = {
            "A": "C=CC(C)C(C)C(C(C)C)C=CC",   # 5-isopropyl-3,4-dimethylocta-1,6-diene
            "B": "C=CCC(C(C)C)C(C)C(C)C=C",   # 4-isopropyl-5,6-dimethylocta-1,7-diene
            "C": "CC=CC(C)C(C(C)C)C=CC",   # 5-isopropyl-3,4-dimethylocta-2,6-diene
            "D": "C=CC(C)C(C)C(C(C)C)CC=C", # 5-isopropyl-3,4-dimethylocta-1,7-diene
        }

        # 3. Define the RCM reaction.
        # This SMARTS pattern finds two alkene groups in a molecule and performs metathesis.
        # It will work for both terminal and internal alkenes.
        rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2].[C:3]=[C:4]>>[C:1]=[C:3].[C:2]=[C:4]')

        correct_option = None

        # 4. Run the simulation for each option.
        for option, smiles in options.items():
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            # Get the two double bonds in the molecule
            double_bonds = [mol.GetBondWithIdx(i) for i in range(mol.GetNumBonds()) if mol.GetBondWithIdx(i).GetBondType() == Chem.BondType.DOUBLE]
            
            # RCM is intramolecular, so we need to apply the reaction to the two double bonds within the same molecule
            if len(double_bonds) == 2:
                # To simulate intramolecular reaction, we break the molecule into two fragments at the double bonds
                # and then run the reaction. This is a standard way to handle this in RDKit.
                ps = rxn.RunReactants((mol,))
                
                # 5. Check the products
                for prods in ps:
                    for p in prods:
                        # We are looking for the cyclic product, not ethene or other small fragments
                        if p.GetNumAtoms() > 4 and Chem.GetSSSR(p) > 0:
                            # Check Constraint 1: Ring size must be 6
                            ring_info = p.GetRingInfo()
                            if ring_info.NumRings() == 1 and len(ring_info.AtomRings()[0]) == 6:
                                # Check Constraint 2: Structure must match target
                                product_canonical_smiles = Chem.MolToSmiles(p, canonical=True)
                                if product_canonical_smiles == target_canonical_smiles:
                                    correct_option = option
                                    break
                    if correct_option:
                        break
            if correct_option:
                break
        
        # 6. Verify the final answer
        given_answer = "D"
        
        if not correct_option:
            return f"Incorrect. The provided answer is {given_answer}, but the code could not find any option that produces the target molecule 5-isopropyl-3,4-dimethylcyclohex-1-ene. There might be an error in the problem statement or options."

        if correct_option == given_answer:
            return "Correct"
        else:
            return f"Incorrect. The provided answer is {given_answer}, but the correct starting material is option {correct_option}. Option {given_answer} does not produce the target molecule 5-isopropyl-3,4-dimethylcyclohex-1-ene."

    except ImportError:
        return "Error: The 'rdkit' library is required to run this check. Please install it using 'pip install rdkit-pypi'."
    except Exception as e:
        return f"An unexpected error occurred: {e}"

# Run the check
result = check_answer()
print(result)