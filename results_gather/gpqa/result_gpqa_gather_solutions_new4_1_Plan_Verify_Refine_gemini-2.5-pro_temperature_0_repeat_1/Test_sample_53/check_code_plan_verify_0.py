import sys
import subprocess

def install_rdkit():
    """
    Installs the RDKit library if it is not already installed.
    RDKit is a cheminformatics toolkit used for molecular manipulation.
    """
    try:
        import rdkit
    except ImportError:
        print("RDKit library not found. Attempting to install...")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
            print("RDKit installed successfully.")
        except Exception as e:
            print(f"Failed to install RDKit: {e}")
            print("Please install it manually using 'pip install rdkit-pypi'")
            sys.exit(1)

def check_chemistry_answer():
    """
    This function verifies the correct starting material for a Ring-Closing Metathesis (RCM) reaction.

    It simulates the RCM for each candidate starting material and compares the resulting
    product with the target molecule to determine the correct answer.
    """
    # Ensure RDKit is available
    install_rdkit()
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # 1. Define the target product and candidate starting materials using SMILES strings.
    
    # Target Product: 5-isopropyl-3,4-dimethylcyclohex-1-ene
    # We generate its canonical SMILES to use as a unique identifier.
    try:
        target_mol = Chem.MolFromSmiles("CC(C)C1CC(C)C(C)C=C1")
        if not target_mol: raise ValueError
        target_canonical_smiles = Chem.MolToSmiles(target_mol, isomericSmiles=True)
    except Exception:
        return "Error: Could not parse the SMILES for the target product '5-isopropyl-3,4-dimethylcyclohex-1-ene'."

    # Candidate starting materials from the multiple-choice options
    options = {
        "A": "CC=C(C)C(C)C(C(C)C)C=CC",  # 5-isopropyl-3,4-dimethylocta-2,6-diene
        "B": "C=CCC(C(C)C)C(C)C(C)C=C",  # 4-isopropyl-5,6-dimethylocta-1,7-diene
        "C": "C=CC(C)C(C)C(C(C)C)C=CC",  # 5-isopropyl-3,4-dimethylocta-1,6-diene
        "D": "C=CC(C)C(C)C(C(C)C)CC=C",  # 5-isopropyl-3,4-dimethylocta-1,7-diene
    }

    # The final answer provided by the LLM to be checked
    llm_answer = "D"

    # 2. Define the RCM reaction using a SMARTS pattern.
    # This pattern finds two double bonds in a molecule and performs metathesis.
    # It creates a new double bond between the internal carbons (forming a ring)
    # and another between the terminal carbons (forming ethene).
    try:
        rxn = AllChem.ReactionFromSmarts('([C:1]=[C:2])~*~([C:3]=[C:4])>>[C:2]=[C:3].[C:1]=[C:4]')
    except Exception:
        return "Error: Could not create the RCM reaction from SMARTS."

    correct_option = None

    # 3. Simulate the RCM reaction for each candidate.
    for option_key, smiles in options.items():
        reactant = Chem.MolFromSmiles(smiles)
        if not reactant:
            continue

        # Run the reaction on the reactant
        products_list = rxn.RunReactants((reactant,))

        # The reaction produces a list of product tuples. We expect one tuple
        # containing two products: the cyclic alkene and ethene.
        if products_list:
            product_tuple = products_list[0]
            for p in product_tuple:
                # We are interested in the major organic product, not the small byproduct (ethene).
                # We can identify the major product by its atom count (ethene has 2 carbons).
                if p.GetNumAtoms() > 2:
                    # Sanitize the molecule and get its canonical SMILES for a reliable comparison.
                    Chem.SanitizeMol(p)
                    product_canonical_smiles = Chem.MolToSmiles(p, isomericSmiles=True)

                    # 4. Compare the simulated product with the target product.
                    if product_canonical_smiles == target_canonical_smiles:
                        correct_option = option_key
                        break  # Found the correct starting material
            if correct_option:
                break # Exit the outer loop as we've found the answer

    # 5. Verify the LLM's answer against the simulation result.
    if correct_option is None:
        return "Failure in verification. The simulation did not find any option that produces the target molecule. There might be an error in the problem statement or the SMILES definitions."

    if correct_option == llm_answer:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer}, but the correct starting material is option {correct_option}. "
                f"The simulation shows that RCM of '{options[correct_option]}' (Option {correct_option}) produces the target molecule, "
                f"while RCM of the other options does not.")

# Execute the verification function and print the result.
result = check_chemistry_answer()
print(result)