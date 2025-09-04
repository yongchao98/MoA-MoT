from rdkit import Chem
from rdkit.Chem import AllChem

def check_organic_synthesis_answer():
    """
    This function programmatically verifies the multi-step organic synthesis problem.
    It uses RDKit to simulate each reaction and compares the final product
    to the proposed answer.
    """
    try:
        # --- 1. Define Molecules and the Proposed Answer ---

        # Starting Material: 3-(hydroxymethyl)-5-(prop-1-en-2-yl)cyclohexan-1-one
        smiles_start = 'O=C1CC(CO)CC(C(C)=C)C1'

        # Options provided in the question
        smiles_options = {
            'A': 'CC(C)C1CC(CO)CCC1',  # (3-isopropylcyclohexyl)methanol
            'B': 'CC(C)C1CC(COCc2ccccc2)CCC1', # (((3-isopropylcyclohexyl)methoxy)methyl)benzene
            'C': 'CC(C)C1CC(CO)CC(C(=NNC(=S(=O)=O)c2ccc(C)cc2)C1)C1', # A plausible structure for the tosylhydrazone intermediate
            'D': 'CCCC1(O)CC(COCc2ccccc2)CC(C(C)C)C1' # 3-((benzyloxy)methyl)-1-butyl-5-isopropylcyclohexan-1-ol
        }
        
        # The final answer to be checked is 'A'
        correct_answer_letter = 'A'

        mol_start = Chem.MolFromSmiles(smiles_start)
        if mol_start is None:
            return "Failed to parse starting material SMILES."

        # --- 2. Simulate the Reaction Sequence Step-by-Step ---

        # Step 1: Williamson Ether Synthesis (Protecting the alcohol)
        # R-OH -> R-OBn
        rxn1 = AllChem.ReactionFromSmarts('[#6]-[O;H1:1]>>[#6]-[O:1]-Cc1ccccc1')
        products1 = rxn1.RunReactants((mol_start,))
        if not products1:
            return "Incorrect. Step 1 (Williamson Ether Synthesis) failed to produce a product."
        mol_step1 = products1[0][0]
        Chem.SanitizeMol(mol_step1)

        # Step 2: Tosylhydrazone Formation (from the ketone)
        # C=O -> C=N-NHTs
        rxn2 = AllChem.ReactionFromSmarts('[C:1]=[O:2]>>[C:1]=NNC(=S(=O)=O)c1ccc(C)cc1')
        products2 = rxn2.RunReactants((mol_step1,))
        if not products2:
            return "Incorrect. Step 2 (Tosylhydrazone Formation) failed to produce a product."
        mol_step2 = products2[0][0]
        Chem.SanitizeMol(mol_step2)

        # Step 3: Shapiro Reaction (Tosylhydrazone to alkene)
        # C(=N-NHTs)-CH2 -> C=C
        rxn3 = AllChem.ReactionFromSmarts('[C:1](=NNC(=S(=O)=O)c1ccc(C)cc1)-[CH2:2]>>[C:1]=[C:2]')
        products3 = rxn3.RunReactants((mol_step2,))
        if not products3:
            return "Incorrect. Step 3 (Shapiro Reaction) failed to produce a product. The SMARTS pattern for the tosylhydrazone reduction did not match."
        mol_step3 = products3[0][0]
        Chem.SanitizeMol(mol_step3)

        # Step 4: Catalytic Hydrogenation and Hydrogenolysis
        mol_step4 = mol_step3
        
        # 4a: Hydrogenate all non-aromatic C=C bonds.
        # We run this in a loop to catch both double bonds.
        rxn4a = AllChem.ReactionFromSmarts('[C;!a:1]=[C;!a:2]>>[C:1]-[C:2]')
        while mol_step4.HasSubstructMatch(Chem.MolFromSmarts('[C;!a]=[C;!a]')):
            products4a = rxn4a.RunReactants((mol_step4,))
            if not products4a: break
            mol_step4 = products4a[0][0]
            Chem.SanitizeMol(mol_step4)

        # Check if the benzyl ether is still present before cleavage
        if not mol_step4.HasSubstructMatch(Chem.MolFromSmiles(smiles_options['B'])):
             return f"Incorrect. After hydrogenation of alkenes (Step 4a), the product should be '{smiles_options['B']}', but it is not. There is an error in the reaction logic."

        # 4b: Cleave the benzyl ether (Hydrogenolysis)
        # R-OBn -> R-OH
        rxn4b = AllChem.ReactionFromSmarts('[O:1]-Cc1ccccc1>>[O:1]H')
        products4b = rxn4b.RunReactants((mol_step4,))
        if not products4b:
            return f"Incorrect. The final step (4b, Hydrogenolysis) failed. The analysis states that the benzyl ether should be cleaved to an alcohol, but this did not happen. This would lead to product B, not A."
        
        mol_final = products4b[0][0]
        Chem.SanitizeMol(mol_final)

        # --- 3. Compare the Final Product with the Proposed Answer ---

        # Get canonical SMILES for comparison (ignores drawing differences)
        final_smiles = Chem.MolToSmiles(mol_final, canonical=True)
        correct_option_mol = Chem.MolFromSmiles(smiles_options[correct_answer_letter])
        correct_option_smiles = Chem.MolToSmiles(correct_option_mol, canonical=True)

        if final_smiles == correct_option_smiles:
            return "Correct"
        else:
            # If the simulation doesn't match the proposed answer, find out what it does match.
            for letter, smiles in smiles_options.items():
                mol = Chem.MolFromSmiles(smiles)
                if Chem.MolToSmiles(mol, canonical=True) == final_smiles:
                    return (f"Incorrect. The provided answer is {correct_answer_letter}, but the simulated reaction "
                            f"sequence leads to the product in option {letter}. The analysis in the provided answer is correct, "
                            f"but the final letter choice is wrong.")
            
            return (f"Incorrect. The provided answer is {correct_answer_letter} ({correct_option_smiles}), but the simulated "
                    f"reaction sequence produces a different molecule with SMILES: {final_smiles}.")

    except ImportError:
        return "RDKit library not found. Please install it using 'pip install rdkit-pypi' to run this check."
    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_organic_synthesis_answer()
print(result)