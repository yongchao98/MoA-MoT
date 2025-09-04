import re

def check_correctness(question, llm_answer):
    """
    Checks the correctness of the LLM's answer to the chemistry question.

    The function encodes the fundamental chemical principles required to solve the problem
    and compares the LLM's reasoning and final answer against this ground truth.
    """

    # 1. Define the ground truth based on chemical principles.
    # Principle 1: Chemoselectivity of the reducing agents.
    reagent_selectivity = {
        'LiBH4': 'ester',  # LiBH4 reduces esters, not acids.
        'BH3': 'acid'      # BH3 reduces acids, not esters.
    }

    # Principle 2: Stereochemical outcome.
    # The reactions do not affect the chiral center, so configuration is retained.
    stereochemistry_rule = 'retained'

    # 2. Deduce the correct answer from the principles.
    # Reaction A: A + LiBH4 -> (R)-product
    # Since stereochemistry is retained, to get an (R) product, A must be (R).
    correct_A_stereo = 'R'

    # Reaction B: B + BH3 -> (S)-product
    # Since stereochemistry is retained, to get an (S) product, B must be (S).
    correct_B_stereo = 'S'

    # Map the deduced configuration to the correct option letter.
    # A) A=(R), B=(R)
    # B) A=(S), B=(R)
    # C) A=(R), B=(S)
    # D) A=(S), B=(S)
    if correct_A_stereo == 'R' and correct_B_stereo == 'S':
        correct_option = 'C'
    else:
        # This part is for robustness, but based on our analysis, it should be 'C'.
        # Logic to determine other options if needed.
        pass

    # 3. Parse the LLM's answer.
    llm_answer_text = llm_answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure: Could not parse the final answer from the LLM's response. Expected format is '<<<X>>>'."
    
    llm_final_choice = match.group(1)

    # 4. Compare the LLM's answer and reasoning to the ground truth.
    # Check 1: Is the final letter correct?
    if llm_final_choice != correct_option:
        return (f"Incorrect. The provided answer is {llm_final_choice}, but the correct answer is {correct_option}.\n"
                f"Reasoning: To get the (R)-product in Reaction A, starting material A must be (R). "
                f"To get the (S)-product in Reaction B, starting material B must be (S). "
                f"This combination (A=R, B=S) corresponds to option C.")

    # Check 2: Is the reasoning correct?
    # We check for keywords that indicate correct understanding.
    reasoning_errors = []
    if not (re.search(r'LiBH₄.*reduce esters', llm_answer_text, re.IGNORECASE)):
        reasoning_errors.append("did not correctly state that LiBH₄ reduces esters.")
    if not (re.search(r'BH₃.*reduce carboxylic acids', llm_answer_text, re.IGNORECASE)):
        reasoning_errors.append("did not correctly state that BH₃ reduces carboxylic acids.")
    if not (re.search(r'retained|conserved|not affect.*chiral center', llm_answer_text, re.IGNORECASE)):
        reasoning_errors.append("did not correctly state that stereochemistry is retained.")
    if not (re.search(r'A must be the \(R\)', llm_answer_text, re.IGNORECASE)):
        reasoning_errors.append("did not correctly deduce that A must be (R).")
    if not (re.search(r'B must be the \(S\)', llm_answer_text, re.IGNORECASE)):
        reasoning_errors.append("did not correctly deduce that B must be (S).")

    if reasoning_errors:
        return f"The final answer {llm_final_choice} is correct, but the reasoning is flawed or incomplete: " + "; ".join(reasoning_errors)

    # If all checks pass, the answer is correct.
    return "Correct"

# The user's provided question and answer to be checked.
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
A chemist performed two reactions by taking two unknown compounds and treated them separately with two different reducing agents. Select the proper starting material for both of the reactions.
A + LiBH4 + H+ ---> (R)-4-ethyltetrahydro-2H-pyran-2-one
B + BH3 + H+ ---> (S)-4-ethyltetrahydro-2H-pyran-2-one

A) A = (R)-3-ethyl-5-isobutoxy-5-oxopentanoic acid , B = (R)-3-ethyl-5-isobutoxy-5-oxopentanoic acid
B) A = (S)-3-ethyl-5-isobutoxy-5-oxopentanoic acid, B = (R)-3-ethyl-5-isobutoxy-5-oxopentanoic acid
C) A = (R)-3-ethyl-5-isobutoxy-5-oxopentanoic acid, B = (S)-3-ethyl-5-isobutoxy-5-oxopentanoic acid
D) A = (S)-3-ethyl-5-isobutoxy-5-oxopentanoic acid , B = (S)-3-ethyl-5-isobutoxy-5-oxopentanoic acid
"""

llm_answer_to_check = """
Here is a step-by-step analysis of the problem to determine the correct answer.

### **Analysis of the Core Chemistry**

The problem asks to identify the correct starting materials, A and B, for two separate reactions that produce enantiomeric products. The key to solving this lies in understanding three main points:
1.  The structure of the starting material and the product.
2.  The chemoselectivity of the two different reducing agents (LiBH₄ and BH₃).
3.  The stereochemical outcome of the reactions.

### **Careful Point 1: Analyzing the Structures**

*   **Starting Material:** The options indicate the starting material is a stereoisomer of `3-ethyl-5-isobutoxy-5-oxopentanoic acid`. This molecule has a five-carbon backbone. At one end (C1) is a carboxylic acid (`-COOH`), and at the other end (C5) is an isobutyl ester (`-COOiBu`). It has a single chiral center at C3, which has an ethyl group attached.
*   **Product:** The product is a stereoisomer of `4-ethyltetrahydro-2H-pyran-2-one`. This is a six-membered cyclic ester (a lactone). It is formed when a `5-hydroxy-3-ethylpentanoic acid` molecule cyclizes.
*   **Correlation:** The chiral center at C3 of the starting material becomes the chiral center at C4 of the lactone product.

### **Careful Point 2: Stereochemistry**

The reactions involve the reduction of a carbonyl group (either the acid or the ester) and a subsequent cyclization. Both of these transformations occur at the ends of the molecule (C1 and C5), which are remote from the chiral center at C3. Since no bonds to the chiral center are made or broken, its absolute configuration (R or S) will be retained throughout the reaction.
*   Therefore, an (R)-configured starting material will produce an (R)-configured lactone.
*   An (S)-configured starting material will produce an (S)-configured lactone.

### **Careful Point 3: Chemoselectivity of the Reducing Agents**

This is the central chemical concept being tested. The starting material has both a carboxylic acid and an ester. The two reducing agents have opposite selectivities for these groups.

**Reaction A: A + LiBH₄ + H⁺ → (R)-lactone**
*   **Reagent:** Lithium borohydride (LiBH₄) is a selective reducing agent. It is known to **reduce esters** to primary alcohols but does not typically reduce carboxylic acids.
*   **Mechanism:** LiBH₄ reduces the isobutyl ester at C5 to a primary alcohol. The carboxylic acid at C1 remains. The resulting intermediate is `5-hydroxy-3-ethylpentanoic acid`. The acidic workup (H⁺) then catalyzes the intramolecular cyclization (lactonization) to form the product.
*   **Conclusion for A:** Since the reaction retains stereochemistry and the product is the **(R)**-lactone, the starting material **A** must be the **(R)**-enantiomer.

**Reaction B: B + BH₃ + H⁺ → (S)-lactone**
*   **Reagent:** Borane (BH₃) is also a selective reducing agent, but with the opposite preference. It is known to **reduce carboxylic acids** to primary alcohols much faster than it reduces esters.
*   **Mechanism:** BH₃ reduces the carboxylic acid at C1 to a primary alcohol. The isobutyl ester at C5 remains. The resulting intermediate is `isobutyl 5-hydroxy-3-ethylpentanoate`. The acidic workup (H⁺) then catalyzes an intramolecular transesterification (the new alcohol attacks the ester) to form the same lactone product.
*   **Conclusion for B:** Since this reaction also retains stereochemistry and the product is the **(S)**-lactone, the starting material **B** must be the **(S)**-enantiomer.

### **Final Conclusion**

Based on the analysis:
*   Starting Material **A** is **(R)-3-ethyl-5-isobutoxy-5-oxopentanoic acid**.
*   Starting Material **B** is **(S)-3-ethyl-5-isobutoxy-5-oxopentanoic acid**.

This combination corresponds to option C.

<<<C>>>
"""

# Run the check
result = check_correctness(question, llm_answer_to_check)
print(result)