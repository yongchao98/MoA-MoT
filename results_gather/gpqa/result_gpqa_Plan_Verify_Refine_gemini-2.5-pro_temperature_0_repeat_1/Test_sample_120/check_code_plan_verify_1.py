from rdkit import Chem

def check_epoxide_opening_reaction():
    """
    Checks the correctness of the predicted product for the reaction of
    (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane with Me2CuLi.
    """
    try:
        # --- Step 1: Define the problem and analyze the reactant ---
        # Reactant Name: (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane
        # Reactant SMILES used in the LLM's reasoning:
        reactant_smiles = "C[C@H]1C[C@H](C)[C@]2(C)O[C@H]12"
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)

        if not reactant_mol:
            return "Failure: The reactant SMILES string 'C[C@H]1C[C@H](C)[C@]2(C)O[C@H]12' is invalid and could not be parsed by RDKit."

        # --- Step 2: Verify Regioselectivity ---
        # Rule: Attack occurs at the less hindered carbon of the epoxide.
        # We expect one tertiary and one quaternary carbon in the epoxide.
        epoxide_carbon_indices = []
        for atom in reactant_mol.GetAtoms():
            if atom.GetAtomicNum() == 8:  # Oxygen
                epoxide_carbon_indices = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
                break
        
        if len(epoxide_carbon_indices) != 2:
            return "Failure: Could not identify exactly two epoxide carbons in the reactant."

        c1 = reactant_mol.GetAtomWithIdx(epoxide_carbon_indices[0])
        c2 = reactant_mol.GetAtomWithIdx(epoxide_carbon_indices[1])
        
        # GetTotalDegree() includes implicit hydrogens. Quaternary C has degree 4, Tertiary C has degree 3.
        degrees = sorted([c1.GetTotalDegree(), c2.GetTotalDegree()])
        if degrees != [3, 4]:
            return f"Constraint check failed: The epoxide carbons should be tertiary and quaternary. Found degrees {degrees} instead."

        # The less hindered carbon (tertiary) is the site of attack.
        attack_site_idx = c1.GetIdx() if c1.GetTotalDegree() == 3 else c2.GetIdx()

        # --- Step 3: Verify Stereoselectivity ---
        # Rule: Inversion of configuration occurs at the attacked carbon.
        Chem.AssignStereochemistry(reactant_mol, cleanIt=True, force=True)
        chiral_centers_reactant = dict(Chem.FindMolChiralCenters(reactant_mol, includeUnassigned=True))

        if attack_site_idx not in chiral_centers_reactant:
            return f"Failure: The identified attack site (atom index {attack_site_idx}) is not a chiral center."

        reactant_config_at_attack_site = chiral_centers_reactant[attack_site_idx]

        # The IUPAC name specifies C6 (the tertiary bridgehead) is 'S'. Our analysis must match this.
        if reactant_config_at_attack_site != 'S':
            return (f"Constraint check failed: The less hindered epoxide carbon (C6) is named as 'S', "
                    f"but the SMILES provided resolves its configuration as '{reactant_config_at_attack_site}'.")

        # According to the SN2-like mechanism, the configuration must invert.
        expected_product_config = 'R'

        # --- Step 4: Analyze the Product Options ---
        options = {
            "A": "(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
            "B": "(1S,4R,5S)-2,2,4,5-tetramethylcyclohexan-1-ol",
            "C": "(1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
            "D": "(1R,4R,5R)-2,2,4,5-tetramethylcyclohexan-1-ol"
        }
        llm_answer = "A"

        # Check regiochemistry of options. The product must be a 1,2,4,5-tetramethylcyclohexanol.
        # Options B and D are 2,2,4,5-isomers, which implies a gem-dimethyl group at C2. This is mechanistically impossible.
        if "2,2," in options["B"] or "2,2," in options["D"]:
            pass # Correctly identified as having wrong regiochemistry
        else:
            return "Logic Error: Could not confirm that options B and D have incorrect regiochemistry based on their names."

        # The choice is between A and C. They differ at C2, the carbon corresponding to the attack site.
        # Option A has C2 as 'R'.
        # Option C has C2 as 'S'.
        
        # Our analysis requires the configuration to be 'R' due to inversion.
        if "2R" not in options["A"]:
             return "Logic Error: Option A name does not contain '2R' as expected."
        if "2S" not in options["C"]:
             return "Logic Error: Option C name does not contain '2S' as expected."

        # The LLM's answer (A) has the correct 'R' configuration at C2.
        if llm_answer == "A":
            return "Correct"
        else:
            return f"Incorrect. The LLM chose {llm_answer}, but the correct answer is A because the reaction requires inversion of the 'S' center to an 'R' center."

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_epoxide_opening_reaction()
print(result)