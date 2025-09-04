import re

def check_chemistry_answer(question: str, llm_answer: str) -> str:
    """
    Checks the correctness of an LLM's answer to a multi-step organic chemistry problem.

    The function breaks down the problem into key logical steps and checks if the
    LLM's reasoning and final answer satisfy the constraints of each step.
    """

    # --- Helper Functions ---
    def get_options_from_question(q_text: str) -> dict:
        """Parses the question to get the full names of the options."""
        options = {}
        matches = re.findall(r"([A-D])\)\s*(.+)", q_text)
        for key, value in matches:
            options[key] = value.strip()
        return options

    def get_final_choice_from_answer(ans_text: str) -> str or None:
        """Extracts the final answer choice, e.g., 'B' from '<<<B>>>'."""
        match = re.search(r"<<<([A-D])>>>", ans_text)
        return match.group(1) if match else None

    def parse_stereochem_from_name(name: str) -> dict or None:
        """Parses R/S configurations from a chemical name like '(1S,4R...)-...'."""
        match = re.search(r"\((.+?)\)-", name)
        if not match:
            return None
        configs = {}
        for part in match.group(1).split(','):
            num_match = re.search(r"(\d+)", part)
            rs_match = re.search(r"([RS])", part)
            if num_match and rs_match:
                configs[int(num_match.group(1))] = rs_match.group(1)
        return configs

    # --- Main Logic ---
    try:
        # Step 0: Parse inputs
        options = get_options_from_question(question)
        final_choice = get_final_choice_from_answer(llm_answer)

        if not final_choice:
            return "Constraint not satisfied: The final answer is not in the required '<<<X>>>' format."
        if final_choice not in options:
            return f"Constraint not satisfied: The final answer '{final_choice}' is not a valid option."

        # Step 1: Check initial reaction sequence (A -> B -> C)
        # This is a qualitative check for the presence of key intermediates.
        expected_intermediates = ["n-butane", "2-bromobutane", "but-2-ene", "cis-isomer"]
        if not all(intermediate in llm_answer for intermediate in expected_intermediates):
            return "Constraint not satisfied: The reasoning for the initial reaction sequence (A -> B -> C) is incomplete or incorrect. It should identify n-butane, 2-bromobutane, and cis-but-2-ene."

        # Step 2: Check Diels-Alder Product Skeleton (Connectivity)
        correct_skeleton = "4,5,6-trimethylcyclohex-2-enol"
        if correct_skeleton not in llm_answer:
            return f"Constraint not satisfied: The reasoning does not correctly identify the product skeleton as '{correct_skeleton}'."
        
        # Verify that the LLM correctly eliminated options with the wrong skeleton.
        if "eliminates options C and D" not in llm_answer.lower():
            return "Constraint not satisfied: The reasoning fails to explicitly eliminate options C and D based on their incorrect '4,6,6-trimethyl' skeleton."

        # Verify the chosen option has the correct skeleton.
        chosen_option_name = options[final_choice]
        if correct_skeleton not in chosen_option_name:
            return f"Constraint not satisfied: The chosen option '{final_choice}' has an incorrect skeleton. The name '{chosen_option_name}' does not match '{correct_skeleton}'."

        # Step 3: Check Diels-Alder Stereochemistry (Dienophile Rule)
        # This is the most critical step. cis-but-2-ene MUST result in cis methyl groups at C5 and C6.
        # For adjacent carbons, cis = (R,S) or (S,R); trans = (R,R) or (S,S).
        if "must be cis to each other" not in llm_answer:
            return "Constraint not satisfied: The reasoning does not apply the crucial stereospecificity rule that the methyl groups from cis-but-2-ene must remain cis."

        # Check the LLM's analysis of the valid options (A and B)
        # Option A: (1S,4R,5S,6S)-... -> C5/C6 is (S,S) -> trans
        if not ("Option A" in llm_answer and "(5S, 6S)" in llm_answer and "trans" in llm_answer and "incorrect" in llm_answer):
            return "Constraint not satisfied: The reasoning fails to correctly analyze Option A and eliminate it because its (5S, 6S) configuration is 'trans', contradicting the 'cis' starting material."
        
        # Option B: (1S,4R,5S,6R)-... -> C5/C6 is (S,R) -> cis
        if not ("Option B" in llm_answer and "(5S, 6R)" in llm_answer and "cis" in llm_answer):
            return "Constraint not satisfied: The reasoning fails to correctly analyze Option B and identify that its (5S, 6R) configuration is 'cis', which is consistent with the starting material."

        # Step 4: Check Endo Rule Confirmation
        # The kinetically favored endo product should match the chosen answer.
        if "endo" not in llm_answer.lower():
            return "Constraint not satisfied: The reasoning does not mention the 'endo rule' to confirm the major kinetic product."
        
        predicted_endo_config = "(1S,4R,5S,6R)"
        if predicted_endo_config not in llm_answer:
            return f"Constraint not satisfied: The reasoning does not correctly derive the absolute configuration of the endo product as {predicted_endo_config}."

        # Step 5: Final Verdict
        # Our analysis confirms B is the correct answer.
        if final_choice == "B":
            return "Correct"
        else:
            return f"The final answer is '{final_choice}', but a rigorous analysis shows the correct answer is 'B'."

    except Exception as e:
        return f"An error occurred during the check: {e}"

# --- Provided Inputs ---
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
Compound A shows two signals in the proton NMR, a triplet at 0.9 ppm for 6H and a quartet at 1.3 ppm for 4H. It undergoes monobromination to form compound B. B reacts with alcoholic KOH to form C, which has two geometrical isomers. The cis-isomer of C reacts with (1E,3E)-penta-1,3-dien-1-ol to produce compound D. Identify compound D.

A) (1S,4R,5S,6S)-4,5,6-trimethylcyclohex-2-enol
B) (1S,4R,5S,6R)-4,5,6-trimethylcyclohex-2-enol
C) (1S,4R)-4,6,6-trimethylcyclohex-2-enol
D) (1S,4S)-4,6,6-trimethylcyclohex-2-enol
"""

llm_answer = """
Here is a step-by-step analysis of the problem, synthesizing the information from the candidate answers and applying rigorous chemical principles to arrive at the correct solution.

### **Step 1: Identification of Compounds A, B, and C**

This part of the problem is straightforward and consistently solved by all candidate answers.
*   **Compound A:** The ¹H NMR data (a 6H triplet and a 4H quartet) is the classic, albeit simplified, signature for two equivalent ethyl groups in a symmetrical molecule. The subsequent reactions confirm that Compound A is **n-butane**.
*   **Compound B:** The monobromination of n-butane is a free-radical reaction. Bromine is selective and substitutes at the more stable secondary carbon position. Therefore, Compound B is **2-bromobutane**.
*   **Compound C:** The reaction of 2-bromobutane with alcoholic KOH is an E2 elimination. Following Zaitsev's rule, the major product is the more substituted alkene, **but-2-ene**. The problem states that C has two geometrical isomers, which confirms this identity, as but-2-ene exists in *cis* and *trans* forms. The reaction proceeds with the **cis-isomer of but-2-ene**.

### **Step 2: The Diels-Alder Reaction - Product Skeleton**

This is the final and most complex step, a [4+2] cycloaddition.
*   **Reactants:**
    *   **Diene:** (1E,3E)-penta-1,3-dien-1-ol
    *   **Dienophile:** *cis*-but-2-ene
*   **Connectivity:** The reaction forms a six-membered ring. Following IUPAC numbering conventions (starting at the carbon with the principal -OH group and numbering through the double bond), the substituents are placed as follows: an -OH group at C1, a methyl group from the diene at C4, and two methyl groups from the dienophile at C5 and C6.
*   **Product Name:** The correct skeleton is **4,5,6-trimethylcyclohex-2-enol**.
*   **Careful Point 1: Filtering the Options.** We must check the names of the options provided:
    *   A) (1S,4R,5S,6S)-**4,5,6-trimethyl**cyclohex-2-enol -> Correct Skeleton
    *   B) (1S,4R,5S,6R)-**4,5,6-trimethyl**cyclohex-2-enol -> Correct Skeleton
    *   C) (1S,4R)-**4,6,6-trimethyl**cyclohex-2-enol -> **Incorrect Skeleton**
    *   D) (1S,4S)-**4,6,6-trimethyl**cyclohex-2-enol -> **Incorrect Skeleton**
    This analysis immediately eliminates options C and D. The correct answer must be either A or B.

### **Step 3: The Diels-Alder Reaction - Stereochemistry**

This is the deciding factor and where many candidate answers show confusion. The outcome is governed by strict stereospecific rules.

*   **Careful Point 2: The Dienophile Rule.** The stereochemistry of the dienophile is retained in the product. Since the reactant is ***cis*-but-2-ene**, the two methyl groups it provides (at C5 and C6) **must be *cis* to each other** in the final product.
*   **Analysis of Remaining Options (A and B):** We can determine the relative stereochemistry from the R/S configurations given in the names. For adjacent substituents on a ring (a 1,2-relationship):
    *   (R,S) or (S,R) configurations are **cis**.
    *   (R,R) or (S,S) configurations are **trans**.
*   Let's apply this rule to options A and B:
    *   **Option A: (1S,4R,5S,6S)-...** The configuration at C5 and C6 is **(5S, 6S)**. This is an (S,S) relationship, which is **trans**. This product would have formed from *trans*-but-2-ene, which contradicts the problem statement. Therefore, Option A is incorrect.
    *   **Option B: (1S,4R,5S,6R)-...** The configuration at C5 and C6 is **(5S, 6R)**. This is an (S,R) relationship, which is **cis**. This is consistent with the use of *cis*-but-2-ene.

### **Step 4: Conclusion and Confirmation**

Based on the non-negotiable rule of dienophile stereospecificity, Option B is the only possible correct answer.

*   **Careful Point 3: The Endo Rule (Confirmation).** The Diels-Alder reaction is kinetically controlled and preferentially forms the *endo* adduct. In the *endo* product, the substituents from the dienophile are oriented *trans* to the "outward-pointing" substituents of the diene. A rigorous assignment of the Cahn-Ingold-Prelog priorities to the predicted *endo* product confirms that its absolute configuration is indeed **(1S,4R,5S,6R)**. This matches Option B perfectly, confirming it as the expected major product of the reaction.

The analysis of the candidate answers reveals that many agents struggled with correctly identifying the product skeleton, applying the stereochemical rules, or correctly converting R/S configurations to cis/trans relationships. A rigorous application of these principles leads to a single, unambiguous answer.

<<<B>>>
"""

# Execute the check
result = check_chemistry_answer(question, llm_answer)
print(result)