import sys
from io import StringIO

# It's good practice to handle potential library import errors.
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    # If rdkit is not installed, the check cannot be performed.
    # We will print an informative message and exit.
    print("Error: RDKit library is not installed. Please install it using 'pip install rdkit-pypi' to run this check.")
    sys.exit(1)

def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer for the chemistry question.
    It determines the structure of the final product and calculates the number of 13C-NMR signals.
    """
    # --- Step 1: Define the final product structure based on the reaction sequence ---
    # The reaction sequence is:
    # 1. Propionaldehyde + EDT -> 2-ethyl-1,3-dithiolane
    # 2. + BuLi -> Lithiated carbanion
    # 3. + Bromoethane -> 2,2-diethyl-1,3-dithiolane
    # 4. + HgCl2/H2O -> 3-pentanone (diethyl ketone)
    # 5. + Wittig from 3-bromopentane -> 3,4-diethylhex-3-ene
    # The SMILES (Simplified Molecular Input Line Entry System) string for the final product,
    # 3,4-diethylhex-3-ene, is CCC(=C(CC)CC)CC.
    final_product_smiles = "CCC(=C(CC)CC)CC"

    # --- Step 2: Use RDKit to analyze the structure ---
    try:
        mol = Chem.MolFromSmiles(final_product_smiles)
        if mol is None:
            return "Error: Could not parse the SMILES string for the final product E."

        # Add hydrogens to get a complete molecular graph for symmetry analysis
        mol_with_hs = Chem.AddHs(mol)

        # --- Step 3: Calculate the number of unique carbon environments ---
        # The number of 13C-NMR signals corresponds to the number of symmetrically
        # non-equivalent carbon atoms. We can find this by ranking the atoms
        # based on their topological environment. Atoms with the same rank are equivalent.
        # The breakTies=False argument ensures that symmetry is fully considered.
        ranks = list(Chem.CanonicalRankAtoms(mol_with_hs, breakTies=False))

        # Filter for carbon atoms (atomic number 6) and find their unique ranks
        carbon_ranks = set()
        for i, atom in enumerate(mol_with_hs.GetAtoms()):
            if atom.GetAtomicNum() == 6:  # It's a Carbon atom
                carbon_ranks.add(ranks[i])

        calculated_signals = len(carbon_ranks)

    except Exception as e:
        return f"An unexpected error occurred during chemical analysis: {e}"

    # --- Step 4: Verify the provided answer against the calculation ---
    # The provided answer's logic correctly concludes there are 3 signals.
    # The question options are: A) 11, B) 6, C) 3, D) 8
    # The provided answer selects option C.

    expected_signals_from_reasoning = 3
    expected_option_for_3_signals = 'C'
    
    # Extract the final choice from the provided answer text
    llm_answer_text = """
    Here is a step-by-step analysis of the problem:

    **Step 1: Analyze the reaction sequence to determine the structure of the final product, E.**

    *   **Reaction 1: Propionaldehyde + EDT / BF₃ → A**
        *   This is a protection reaction. The aldehyde group of propionaldehyde (CH₃CH₂CHO) reacts with 1,2-ethanedithiol (EDT) to form a cyclic thioacetal, also known as a 1,3-dithiolane.
        *   Product **A** is 2-ethyl-1,3-dithiolane.

    *   **Reaction 2: A + BuLi → B**
        *   This is an acid-base reaction. n-Butyllithium (BuLi) is a very strong base that removes the most acidic proton from A. The proton on the carbon between the two sulfur atoms is acidic (pKa ≈ 31) because the resulting negative charge is stabilized by the adjacent sulfur atoms. This is an example of umpolung (polarity inversion).
        *   Product **B** is the nucleophilic lithium salt of the 2-ethyl-1,3-dithiolane carbanion.

    *   **Reaction 3: B + Bromoethane → C**
        *   This is a nucleophilic substitution (Sₙ2) reaction. The carbanion (B) attacks the electrophilic carbon of bromoethane (CH₃CH₂Br), forming a new C-C bond.
        *   Product **C** is 2,2-diethyl-1,3-dithiolane.

    *   **Reaction 4: C + HgCl₂ / H₂O / H⁺ → D**
        *   This is a deprotection reaction. The thioacetal is hydrolyzed back to a carbonyl group using mercuric chloride and aqueous acid.
        *   Product **D** is 3-pentanone (also known as diethyl ketone), (CH₃CH₂)₂C=O.

    *   **Reaction 5: D + PPh₃ / 3-bromopentane / BuLi → E**
        *   This is a Wittig reaction, which converts a ketone into an alkene.
        *   First, the Wittig reagent (ylide) is formed from 3-bromopentane ((CH₃CH₂)₂CHBr), triphenylphosphine (PPh₃), and BuLi. The ylide is (CH₃CH₂)₂C=PPh₃.
        *   Then, the ylide reacts with the ketone (D, 3-pentanone). The oxygen atom of the ketone is replaced by the carbon group from the ylide.
        *   (CH₃CH₂)₂C=O + (CH₃CH₂)₂C=PPh₃ → (CH₃CH₂)₂C=C(CH₂CH₃)₂ + Ph₃P=O
        *   The final product, **E**, is 3,4-diethylhex-3-ene.

    **Step 2: Analyze the structure of E to determine the number of ¹³C-NMR signals.**

    *   **Structure of E:** 3,4-diethylhex-3-ene
        ```
              CH₂CH₃
              |
        CH₃CH₂-C=C-CH₂CH₃
              |
              CH₂CH₃
        ```
    *   **Symmetry Analysis:** The molecule is highly symmetrical. It has a C₂ axis of rotation perpendicular to the C=C double bond and a plane of symmetry that bisects the double bond.
    *   **Identifying Unique Carbons:** Due to this symmetry, several groups of carbons are chemically equivalent.
        1.  **Alkene Carbons (C=C):** The two carbons of the double bond are equivalent to each other. This will produce **one signal**.
        2.  **Methylene Carbons (-CH₂-):** All four ethyl groups are equivalent. Therefore, the four methylene carbons are also equivalent to each other. This will produce a **second signal**.
        3.  **Methyl Carbons (-CH₃):** Similarly, the four terminal methyl carbons of the ethyl groups are all equivalent. This will produce a **third signal**.

    **Conclusion:**

    The final product, 3,4-diethylhex-3-ene, has three distinct types of carbon atoms. Therefore, its ¹³C-NMR spectrum will show **3** signals. This corresponds to option C.

    <<<C>>>
    """
    final_choice = llm_answer_text.strip().split('<<<')[-1].split('>>>')[0]
    
    # Check 1: Does the calculated number of signals match the number from the reasoning?
    if calculated_signals != expected_signals_from_reasoning:
        return f"Incorrect. The provided answer's reasoning states there are {expected_signals_from_reasoning} signals, but a computational check on the final product (3,4-diethylhex-3-ene) reveals there are {calculated_signals} signals."

    # Check 2: Does the final chosen option match the correct option based on the reasoning?
    if final_choice != expected_option_for_3_signals:
        return f"Incorrect. The reasoning correctly identifies {expected_signals_from_reasoning} signals, which corresponds to option '{expected_option_for_3_signals}'. However, the final answer provided is '<<<{final_choice}>>>'."

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)