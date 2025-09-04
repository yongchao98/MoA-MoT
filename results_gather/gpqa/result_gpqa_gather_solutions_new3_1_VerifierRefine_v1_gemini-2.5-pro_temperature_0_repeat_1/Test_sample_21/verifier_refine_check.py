import sys
import io

def check_correctness_of_electrochemistry_answer():
    """
    This function checks the correctness of the provided answer to the electrochemistry question.
    It verifies the two key scientific facts required to answer the question.
    """
    try:
        # --- Ground Truth Data from Electrochemistry ---

        # 1. Thermodynamic Aspect: Oxidizing strength is determined by standard reduction potential (E°).
        # A higher E° means a stronger oxidant.
        # We compare the potential in basic solution to the potential in acidic solution.
        E_potential_acidic = 1.23  # Volts for O₂(g) + 4H⁺(aq) + 4e⁻ → 2H₂O(l)
        E_potential_basic = 0.40   # Volts for O₂(g) + 2H₂O(l) + 4e⁻ → 4OH⁻(aq)

        # 2. Kinetic Aspect: Reaction rate.
        # The oxygen reduction reaction (ORR) is known to be kinetically sluggish.
        # When comparing media, the reaction is generally slower in acidic solutions than in basic solutions.
        # The question asks about the rate in acidic solutions.
        kinetic_rate_in_acidic = "slower"

        # --- Provided Answer and Options ---
        options = {
            "A": "stronger - slower",
            "B": "stronger - faster",
            "C": "weaker - faster",
            "D": "weaker - slower"
        }
        
        # The final answer provided by the LLM being checked.
        llm_answer_key = "D"
        llm_answer_text = options.get(llm_answer_key)

        if not llm_answer_text:
            return f"Invalid answer key '{llm_answer_key}'. The key must be one of {list(options.keys())}."

        # --- Verification Logic ---

        # Step 1: Verify the thermodynamic conclusion.
        # "Thermodynamically, oxygen is a ...... oxidant in basic solutions."
        if E_potential_basic < E_potential_acidic:
            correct_thermo_term = "weaker"
        else:
            correct_thermo_term = "stronger"

        # Step 2: Verify the kinetic conclusion.
        # "Kinetically, oxygen reacts ...... in acidic solutions."
        correct_kinetic_term = kinetic_rate_in_acidic

        # Step 3: Construct the correct full answer string.
        correct_combination = f"{correct_thermo_term} - {correct_kinetic_term}"

        # Step 4: Compare the LLM's answer with the correct combination.
        llm_thermo_term, llm_kinetic_term = llm_answer_text.split(" - ")

        error_messages = []
        
        # Check thermodynamic part
        if llm_thermo_term != correct_thermo_term:
            error_messages.append(
                f"Constraint 1 (Thermodynamics) is not satisfied. "
                f"The standard reduction potential of oxygen in basic solution (E°=+{E_potential_basic}V) is lower than in acidic solution (E°=+{E_potential_acidic}V). "
                f"This means oxygen is a '{correct_thermo_term}' oxidant in basic solutions, but the answer claims it is '{llm_thermo_term}'."
            )

        # Check kinetic part
        if llm_kinetic_term != correct_kinetic_term:
            error_messages.append(
                f"Constraint 2 (Kinetics) is not satisfied. "
                f"The oxygen reduction reaction is known to be kinetically '{correct_kinetic_term}' in acidic solutions compared to basic solutions. "
                f"The answer incorrectly claims it is '{llm_kinetic_term}'."
            )

        if error_messages:
            return "\n".join(error_messages)
        else:
            return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_correctness_of_electrochemistry_answer()
print(result)