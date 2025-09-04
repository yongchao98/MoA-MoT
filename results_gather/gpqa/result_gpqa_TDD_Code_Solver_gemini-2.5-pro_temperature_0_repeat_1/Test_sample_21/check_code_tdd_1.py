import re

def check_electrochem_answer(llm_answer: str):
    """
    Checks the correctness of the LLM's answer for the given electrochemistry question.

    The question is:
    A student regrets that he fell asleep during a lecture in electrochemistry, facing the following incomplete statement in a test:
    Thermodynamically, oxygen is a …… oxidant in basic solutions. Kinetically, oxygen reacts …… in acidic solutions.
    Which combination of weaker/stronger and faster/slower is correct?

    A) stronger – faster
    B) weaker – faster
    C) weaker - slower
    D) stronger – slower
    """
    
    # Step 1: Determine the correct answer based on scientific principles.
    # Thermodynamics: The standard reduction potential of oxygen is E° = +0.40 V in basic solution
    # and E° = +1.23 V in acidic solution. Since 0.40 V < 1.23 V, oxygen is a WEAKER oxidant in basic solutions.
    # Kinetics: The reduction of O₂ has a high activation energy due to the strong O=O bond and the
    # multi-electron transfer mechanism, making the reaction kinetically SLOW.
    # Therefore, the correct combination is "weaker - slower".
    correct_option = 'C'
    possible_options = {'A', 'B', 'C', 'D'}

    # Step 2: Parse the LLM's answer to find the selected option.
    # The provided answer is a conversational text that does not select an option.
    # We will search for an option letter to confirm its absence.
    # We use word boundaries (\b) to find standalone letters A, B, C, or D.
    
    found_options = re.findall(r'\b([A-D])\b', llm_answer.upper())
    
    # Also check for the specific format <<<X>>>
    match = re.search(r'<<<([A-D])>>>', llm_answer.upper())
    if match:
        found_options.append(match.group(1))

    unique_found_options = set(found_options)

    if not unique_found_options:
        return "The answer is incorrect because it does not address the question. The response is a conversational statement and fails to select any of the multiple-choice options (A, B, C, or D)."
    
    if len(unique_found_options) > 1:
        return f"The answer is incorrect because it is ambiguous, mentioning multiple options: {', '.join(sorted(list(unique_found_options)))}."

    extracted_answer = unique_found_options.pop()

    # Step 3: Compare the extracted answer with the correct answer.
    if extracted_answer == correct_option:
        return "Correct"
    else:
        return (f"The answer is incorrect. It selected option {extracted_answer}, but the correct option is {correct_option}. "
                f"Reasoning: Oxygen is a 'weaker' oxidant in basic solution (E° = +0.40 V) compared to acidic solution (E° = +1.23 V), "
                f"and its reduction is kinetically 'slower' due to a high activation energy.")

# The provided LLM answer to check:
llm_response_to_check = "Excellent, the previous answer was correct. I'm ready for the next question. Please provide it, and I will follow the Test-Driven Development process to find the solution."

# Execute the check and print the result.
result = check_electrochem_answer(llm_response_to_check)
print(result)