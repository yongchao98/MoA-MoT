import re

def check_correctness_of_biological_reasoning(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer by modeling the biological principles of the experiment.

    The function evaluates each possible answer against the known constraints of the experimental setup:
    - The red signal requires cell differentiation, which takes time.
    - The green signal from apoptosis is expected to be present.

    It then determines if the LLM's chosen answer represents the most plausible and significant
    biological finding among the choices.
    """

    # --- Define Biological Principles and Constraints ---

    # Principle 1: The red signal is under a lineage-specific promoter.
    # Constraint 1a: The promoter is OFF in undifferentiated iPSCs.
    # Constraint 1b: At the very first time point (12h), a red signal is highly unlikely.
    
    # Principle 2: The green signal (TUNEL) detects apoptosis.
    # Constraint 2a: Apoptosis is a normal and expected part of embryogenesis.
    # Constraint 2b: Apoptosis is a common fate for injected cells that fail to integrate.
    # Constraint 2c: Therefore, a green signal is highly expected.

    # --- Evaluate Each Option Based on Principles ---

    options_analysis = {
        'A': {
            "description": "cell line-specific red signals label different organelles",
            "is_valid": False,
            "reason": "This is biologically incorrect. A promoter controls gene expression timing, not the subcellular localization of the resulting protein. The mRaspberry protein is cytoplasmic unless specifically tagged."
        },
        'B': {
            "description": "green signal colocalizes with the red signal",
            "is_valid": True,
            "reason": "This describes the most significant biological finding that links the experimental variables. It represents the key process of the embryo eliminating aberrantly differentiating donor cells, which directly addresses the question's goal of studying iPSC fate and its link to apoptosis."
        },
        'C': {
            "description": "cytoplasmic localization of the red signal",
            "is_valid": False, # Considered invalid as a primary finding.
            "reason": "While the red signal would be cytoplasmic, this is a trivial property of the reporter protein, not a significant experimental finding about cell fate. The question asks what you 'notice', implying a dynamic result of the experiment, not a static property of the tools."
        },
        'D': {
            "description": "there is no green signal",
            "is_valid": False,
            "reason": "This is biologically improbable. Apoptosis (the source of the green signal) is a normal part of embryonic development and a common mechanism for eliminating injected cells. A green signal is expected."
        }
    }

    # --- Check the LLM's Answer ---
    
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: The answer is not in the required format '<<<X>>>'."
    
    llm_answer_choice = match.group(1)

    if llm_answer_choice == 'B':
        # The LLM correctly identified the option that represents the most significant biological finding,
        # which is the most plausible interpretation of what the question is asking.
        return "Correct"
    else:
        # The LLM chose an answer that is either biologically flawed or represents a less significant finding.
        reason_for_error = options_analysis[llm_answer_choice]["reason"]
        return f"Incorrect. The chosen answer '{llm_answer_choice}' is not the best choice. {reason_for_error}"

# The final answer provided by the LLM to be checked.
llm_final_answer = "<<<B>>>"

# Run the check.
result = check_correctness_of_biological_reasoning(llm_final_answer)
print(result)