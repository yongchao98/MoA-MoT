import re

def check_organic_chemistry_answer(llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given organic chemistry question.

    The question is:
    What is the major outcome of the reaction between 4,4-dimethylcyclopent-1-enol and bromine?

    A) 4-bromo-4,4-dimethylcyclopentanone
    B) (1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol
    C) 2-bromo-4,4-dimethylcyclopentanone
    D) (1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol
    """

    # Step 1: Define the correct answer based on chemical principles.
    # The reaction is between an enol (4,4-dimethylcyclopent-1-enol) and a halogen (Br2).
    # This is a classic alpha-halogenation reaction.
    # The enol's double bond (at C1-C2) is nucleophilic, especially at the alpha-carbon (C2).
    # The alpha-carbon attacks the electrophilic bromine.
    # The molecule then tautomerizes to the more stable keto form.
    # This means the hydroxyl at C1 becomes a ketone (C=O), and a bromine atom is added to C2.
    # The resulting product is 2-bromo-4,4-dimethylcyclopentanone.
    
    options = {
        "A": "4-bromo-4,4-dimethylcyclopentanone",
        "B": "(1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol",
        "C": "2-bromo-4,4-dimethylcyclopentanone",
        "D": "(1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol"
    }
    
    correct_option = "C"
    correct_product_name = options[correct_option]

    # Step 2: Extract the final answer from the LLM's response.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    
    if not match:
        return "The answer is malformed. It does not contain the final answer in the format <<<X>>>."

    submitted_option = match.group(1)

    # Step 3: Compare the submitted answer with the correct answer and provide a reason if incorrect.
    if submitted_option == correct_option:
        return "Correct"
    else:
        submitted_product_name = options.get(submitted_option, "an invalid option")
        reason = (
            f"The provided answer is incorrect. It selected option {submitted_option} ('{submitted_product_name}').\n"
            f"The correct answer is option {correct_option} ('{correct_product_name}').\n\n"
            f"**Reasoning:** The reaction of an enol (4,4-dimethylcyclopent-1-enol) with bromine (Br2) is an alpha-halogenation. "
            f"This reaction is driven by the formation of the thermodynamically stable ketone product. "
            f"The bromine atom adds to the alpha-carbon (C2), and the enol group at C1 tautomerizes to a ketone. "
            f"This results in the formation of 2-bromo-4,4-dimethylcyclopentanone, which corresponds to option C."
        )
        return reason

# The full response from the LLM to be checked
llm_full_response = """
Here is a step-by-step analysis to determine the final answer based on the provided candidate responses.

**1. Analysis of the Chemical Reaction**

*   **Reactants**: The reaction involves 4,4-dimethylcyclopent-1-enol and bromine (Br₂).
*   **Key Functional Group**: The starting organic molecule is an **enol**. An enol features a hydroxyl group (-OH) attached directly to a carbon atom within a carbon-carbon double bond (C=C-OH).
*   **Reaction Type**: The reaction of an enol with a halogen like bromine is a fundamental and well-established reaction in organic chemistry.
*   **Dominant Pathway**: There are two main possibilities: simple addition across the double bond or alpha-substitution. The vast majority of the candidate answers correctly identify that the major pathway is **alpha-substitution** (also known as alpha-halogenation).
*   **Reasoning (Thermodynamic Driving Force)**: This pathway is strongly favored because the final product contains a carbonyl group (C=O). The carbon-oxygen double bond is significantly more stable than the carbon-carbon double bond of the starting enol, providing a powerful thermodynamic driving force for the reaction to proceed in this direction.

**2. Prediction of the Product Structure**

The mechanism of alpha-substitution proceeds as follows:
*   The enol's double bond acts as a nucleophile, attacking the electrophilic bromine.
*   The bromine atom attaches to the "alpha-carbon" – the carbon of the double bond that is *not* attached to the hydroxyl group (C2 in this case).
*   The molecule then rapidly rearranges (tautomerizes) to its more stable keto form. The hydroxyl group at C1 becomes a ketone (C=O).
*   The dimethyl groups at C4 are spectators and remain unchanged.

This process results in the formation of **2-bromo-4,4-dimethylcyclopentanone**.

**3. Evaluation of the Candidate Answers**

*   **Consensus on Reasoning**: There is an overwhelming consensus among the candidate answers (at least 15 out of 17) regarding the chemical reasoning. They correctly identify the reaction as an alpha-halogenation of an enol and the product as 2-bromo-4,4-dimethylcyclopentanone.
*   **Inconsistency in Final Choice**: The final letter choices (<<<A>>>, <<<B>>>, <<<C>>>, <<<D>>>) are highly inconsistent across the different answers. A close look reveals that the options presented within each candidate's response were scrambled. For example, in some answers, "2-bromo-4,4-dimethylcyclopentanone" was option A, in others B, and in others C. This means the final letter from each candidate is unreliable.
*   **Conclusion from Analysis**: The most reliable information is the consistent, correct chemical reasoning shared by the vast majority of the agents. This reasoning points unequivocally to a single product.

**4. Final Determination**

Based on the solid consensus of chemical principles, the major product of the reaction is **2-bromo-4,4-dimethylcyclopentanone**. Matching this name to the options provided in the question:

*   A) 4-bromo-4,4-dimethylcyclopentanone
*   B) (1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol
*   **C) 2-bromo-4,4-dimethylcyclopentanone**
*   D) (1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol

The correct option is C.

<<<C>>>
"""

# Execute the check
result = check_organic_chemistry_answer(llm_full_response)
print(result)