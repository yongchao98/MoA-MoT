import re

def check_chemistry_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the selected option for a series of Michael addition reactions.

    The function validates the proposed products (A, B) and reactant (C) against
    the known outcomes of the specified chemical reactions.

    Args:
        llm_answer_text: A string containing the LLM's answer, expected to be
                         in the format <<<X>>> where X is one of A, B, C, or D.

    Returns:
        A string indicating "Correct" if the answer is right, or a detailed
        explanation of the error(s) if it is incorrect.
    """

    # --- Define Correct Chemical Outcomes ---

    def get_correct_A():
        """Determines the correct product for Reaction A."""
        return {
            "name": "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate",
            "explanation": "In Reaction A, the enolate of dimethyl malonate attacks the β-carbon of methyl (E)-3-(p-tolyl)acrylate. The resulting product, after protonation, is correctly named trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate based on its structure: (p-tolyl)-CH(CH(COOCH₃)₂) - CH₂-COOCH₃."
        }

    def get_correct_B():
        """Determines the correct product for Reaction B."""
        return {
            "name": "3-(2-oxocyclohexyl)butanenitrile",
            "explanation": "Reaction B is a Stork enamine synthesis. After the Michael addition, the acidic workup (H₃O⁺) hydrolyzes the intermediate iminium salt back to a ketone. The major product is the thermodynamically more stable keto tautomer, 3-(2-oxocyclohexyl)butanenitrile, not the enol form."
        }

    def get_correct_C():
        """Determines the correct reactant for Reaction C."""
        return {
            "name": "cyclohexane-1,3-dione",
            "explanation": "Reaction C is a Michael addition where the product is given. A retrosynthetic analysis shows that the (3-oxobutyl) group comes from but-3-en-2-one. This group is attached to the C2 position of the other reactant. The highly acidic C2 position of cyclohexane-1,3-dione makes it the correct Michael donor (reactant C). The name '2-hydroxycyclohexane-1,3-dione' is not a standard representation for this reactant."
        }

    # --- Parse and Validate the LLM's Answer ---

    # The options as defined in the question
    options = {
        'A': {
            'A': 'trimethyl 3-(p-tolyl)propane-1,1,2-tricarboxylate',
            'B': '3-(2-hydroxycyclohex-1-en-1-yl)butanenitrile',
            'C': 'cyclohexane-1,3-dione'
        },
        'B': {
            'A': 'trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate',
            'B': '3-(2-oxocyclohexyl)butanenitrile',
            'C': '2-hydroxycyclohexane-1,3-dione'
        },
        'C': {
            'A': 'trimethyl 3-(p-tolyl)propane-1,1,2-tricarboxylate',
            'B': '3-(2-hydroxycyclohex-1-en-1-yl)butanenitrile',
            'C': '2-hydroxycyclohexane-1,3-dione'
        },
        'D': {
            'A': 'trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate',
            'B': '3-(2-oxocyclohexyl)butanenitrile',
            'C': 'cyclohexane-1,3-dione'
        }
    }

    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. The answer should be in the format <<<X>>> where X is A, B, C, or D."

    candidate_option_key = match.group(1)
    proposed_answer = options.get(candidate_option_key)

    if not proposed_answer:
        # This case should not be reached if the regex is correct, but it's good practice.
        return f"Invalid option '{candidate_option_key}' provided."

    # Get the scientifically correct answers
    correct_A = get_correct_A()
    correct_B = get_correct_B()
    correct_C = get_correct_C()

    errors = []

    # Check component A
    if proposed_answer['A'] != correct_A['name']:
        errors.append(f"Component A is incorrect. The proposed product is '{proposed_answer['A']}', but the correct product is '{correct_A['name']}'.\nReason: {correct_A['explanation']}")

    # Check component B
    if proposed_answer['B'] != correct_B['name']:
        errors.append(f"Component B is incorrect. The proposed product is '{proposed_answer['B']}', but the correct product is '{correct_B['name']}'.\nReason: {correct_B['explanation']}")

    # Check component C
    if proposed_answer['C'] != correct_C['name']:
        errors.append(f"Component C is incorrect. The proposed reactant is '{proposed_answer['C']}', but the correct reactant is '{correct_C['name']}'.\nReason: {correct_C['explanation']}")

    if not errors:
        return "Correct"
    else:
        return "\n\n".join(errors)

# The final answer provided by the LLM to be checked
llm_final_answer = """
Here is a step-by-step analysis of the three reactions to determine the correct option.

### **Reaction A: dimethyl malonate + methyl (E)-3-(p-tolyl)acrylate + (NaOEt, EtOH) --> (A)**

1.  **Identify the Nucleophile and Electrophile:**
    *   **Nucleophile (Michael Donor):** Dimethyl malonate is a classic Michael donor. The strong base, sodium ethoxide (NaOEt), deprotonates the acidic α-carbon between the two ester groups, creating a resonance-stabilized enolate, `⁻CH(COOCH₃)₂`.
    *   **Electrophile (Michael Acceptor):** Methyl (E)-3-(p-tolyl)acrylate is an α,β-unsaturated ester. The β-carbon (the one attached to the p-tolyl group) is electrophilic.

2.  **Determine the Product Structure:**
    *   The malonate enolate attacks the β-carbon of the acrylate in a conjugate addition.
    *   The resulting enolate intermediate is protonated by the ethanol solvent.
    *   The final product (A) has the structure: `(p-tolyl)-CH(CH(COOCH₃)₂) - CH₂-COOCH₃`.

3.  **Name the Product (A):**
    *   To match the options, we use the `propane-tricarboxylate` naming system. The main three-carbon chain is the one connecting the carboxylate groups.
    *   The structure can be viewed as `(CH₃OOC)₂CH(1)-CH(2)(p-tolyl)-CH₂(3)-COOCH₃`.
    *   This is a propane chain with two methoxycarbonyl groups at position 1, a p-tolyl group at position 2, and the third methoxycarbonyl group at position 3.
    *   The correct name is **trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate**. This eliminates options A and C.

### **Reaction B: 1-(cyclohex-1-en-1-yl)piperidine + (E)-but-2-enenitrile + (MeOH, H₃O⁺) --> (B)**

1.  **Identify the Nucleophile and Electrophile:**
    *   **Nucleophile (Michael Donor):** 1-(cyclohex-1-en-1-yl)piperidine is an enamine derived from cyclohexanone. In this Stork enamine synthesis, the enamine acts as a carbon nucleophile at the α-position of the original ketone.
    *   **Electrophile (Michael Acceptor):** (E)-but-2-enenitrile (crotononitrile) is an α,β-unsaturated nitrile.

2.  **Determine the Product Structure:**
    *   The enamine attacks the β-carbon of the nitrile.
    *   The subsequent workup with aqueous acid (H₃O⁺) hydrolyzes the intermediate iminium salt back to a ketone.
    *   The final product (B) is a cyclohexanone ring substituted at the α-carbon with a `-CH(CH₃)CH₂CN` group.

3.  **Name the Product (B):**
    *   The major product is the thermodynamically stable keto tautomer, not the enol form. Therefore, the name will contain "oxocyclohexyl".
    *   Naming the molecule as a derivative of butanenitrile (the principal functional group), the substituent is a `(2-oxocyclohexyl)` group attached to position 3 of the butanenitrile chain.
    *   The correct name is **3-(2-oxocyclohexyl)butanenitrile**. This is consistent with options B and D.

### **Reaction C: C + but-3-en-2-one + (KOH, H₂O) ---> 2-(3-oxobutyl)cyclohexane-1,3-dione**

1.  **Identify the Reactant C via Retrosynthesis:**
    *   **Product:** 2-(3-oxobutyl)cyclohexane-1,3-dione.
    *   **Michael Acceptor:** but-3-en-2-one (methyl vinyl ketone, MVK).
    *   **Analysis:** The `(3-oxobutyl)` group, `-CH₂CH₂COCH₃`, in the product comes from the MVK acceptor. This group is attached to the C2 position of the other reactant. The C2 position of a 1,3-dione is highly acidic and is the site of deprotonation by the base (KOH) to form the nucleophilic enolate.
    *   Therefore, the starting Michael donor (C) must be **cyclohexane-1,3-dione**.

### **Conclusion**

Based on the step-by-step analysis:
*   **A** = trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate
*   **B** = 3-(2-oxocyclohexyl)butanenitrile
*   **C** = cyclohexane-1,3-dione

This combination of structures correctly matches all parts of option B (and the identical option D).

<<<B>>>
"""

# Run the check
result = check_chemistry_answer(llm_final_answer)
print(result)