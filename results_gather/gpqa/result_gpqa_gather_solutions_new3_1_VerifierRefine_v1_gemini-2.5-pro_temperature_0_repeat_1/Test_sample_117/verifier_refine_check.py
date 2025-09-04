import re

def check_correctness(question, final_answer):
    """
    Checks the correctness of the final answer for the given organic chemistry question.

    Args:
        question (str): The question text (not used in this specific implementation but good practice).
        final_answer (str): The final answer provided by the LLM, in the format "<<<X>>>".

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """
    # Extract the letter from the final answer format, e.g., "<<<C>>>" -> "C"
    match = re.search(r'<<<([A-D])>>>', final_answer)
    if not match:
        return f"Invalid answer format. Expected format like '<<<A>>>', but got '{final_answer}'."
    
    selected_option = match.group(1)

    # Define the correct answer and the reasoning based on organic chemistry principles.
    # Reaction: 4,4-dimethylcyclopent-1-enol + Br2
    # This is a reaction of an enol with a halogen.
    # The major pathway is alpha-halogenation, not electrophilic addition.
    # The enol's double bond (C1=C2) attacks Br2. Bromine adds to the alpha-carbon (C2).
    # The enol tautomerizes to the more stable keto form (C=O at C1).
    # Product: 2-bromo-4,4-dimethylcyclopentanone.
    
    correct_option = 'C'
    
    # Check if the selected option is the correct one.
    if selected_option == correct_option:
        return "Correct"
    else:
        # Provide specific reasons why other options are incorrect.
        if selected_option in ['A', 'B']:
            return (f"Incorrect. The selected answer corresponds to a 1,2-dibromo-4,4-dimethylcyclopentanol. "
                    "This is a product of simple electrophilic addition. While possible for a standard alkene, "
                    "the reaction of an enol with a halogen is dominated by the alpha-substitution pathway, "
                    "which is driven by the formation of the much more stable carbonyl (C=O) group. "
                    "The major product is an alpha-bromo ketone, not a dibromo-alcohol.")
        elif selected_option == 'D':
            return (f"Incorrect. The selected answer is 4-bromo-4,4-dimethylcyclopentanone. "
                    "The reaction occurs at the alpha-position (C2) relative to the carbonyl group that is formed. "
                    "Position C4 is a quaternary carbon (it has no hydrogen atoms) and is not the electronically activated site for this reaction.")
        else:
            # This case should not be hit with the given options, but it's good for robustness.
            return f"An invalid option '{selected_option}' was chosen."

# The user's provided question and the final answer to be checked.
question_text = "What is the major outcome of the reaction between 4,4-dimethylcyclopent-1-enol and bromine?"
final_answer_from_llm = "<<<C>>>"

# Run the check.
result = check_correctness(question_text, final_answer_from_llm)
print(result)