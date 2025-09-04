import re

def check_correctness(llm_answer):
    """
    Checks the correctness of the LLM's answer for the impulse approximation question.

    The user's question asks which set of assumptions jointly imply the "impulse" approximation in nuclear physics.
    The correct assumptions are:
    1. The interaction is with a single nucleon (one-body current).
    2. The struck nucleon leaves the nucleus without further interaction (no final-state interactions).
    4. The interaction is so rapid that the nucleon's binding forces are negligible during the event.

    Assumption 3, regarding the internal quark structure of the nucleon, is a separate physical consideration and not part of the standard impulse approximation.

    Therefore, the correct combination is (1, 2, 4), which corresponds to option A.
    """
    try:
        # Use regex to find the answer in the format <<<ANSWER>>>
        match = re.search(r'<<<(.+?)>>>', llm_answer)
        if not match:
            return "The answer is not provided in the required format '<<<...>>>'."

        # Extract and clean the answer
        answer = match.group(1).strip().upper()

        # Define the correct option
        correct_option = 'A'

        # Check if the extracted answer is correct
        if answer == correct_option:
            return "Correct"
        else:
            # Provide specific feedback for incorrect answers
            if answer == 'B':
                return "Incorrect. Option B includes assumption 3 (non-relativistic quarks). The impulse approximation concerns the nucleon's behavior within the nucleus, not its internal quark structure."
            elif answer == 'C':
                return "Incorrect. Option C includes assumption 3 (non-relativistic quarks). The impulse approximation concerns the nucleon's behavior within the nucleus, not its internal quark structure."
            elif answer == 'D':
                return "Incorrect. Option D omits assumption 1 (interaction with a single nucleon), which is a fundamental part of the impulse approximation. It also incorrectly includes assumption 3."
            else:
                return f"Incorrect. The provided answer '{answer}' is not a valid option. The valid options are A, B, C, or D."

    except Exception as e:
        return f"An error occurred during checking: {e}"
