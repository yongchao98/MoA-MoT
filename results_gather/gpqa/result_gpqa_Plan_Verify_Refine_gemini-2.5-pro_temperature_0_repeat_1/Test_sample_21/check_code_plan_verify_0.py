import re

def check_answer_correctness():
    """
    Checks the correctness of the answer to the electrochemistry question.

    The function verifies two key concepts:
    1. Thermodynamic strength of O2 as an oxidant in basic vs. acidic solution.
    2. Kinetic rate of the O2 reduction reaction.
    """

    # --- Part 1: Thermodynamic Analysis ---
    # Standard reduction potentials (E°) determine thermodynamic strength.
    # A higher E° means a stronger oxidant.
    E_standard_acidic = 1.23  # Volts for O2 + 4H+ + 4e- -> 2H2O
    E_standard_basic = 0.40   # Volts for O2 + 2H2O + 4e- -> 4OH-

    # The question asks about O2 in basic solution. We compare its potential to the acidic case.
    if E_standard_basic < E_standard_acidic:
        thermodynamic_conclusion = "weaker"
    else:
        thermodynamic_conclusion = "stronger"

    # --- Part 2: Kinetic Analysis ---
    # The Oxygen Reduction Reaction (ORR) is known to be kinetically slow due to a high activation energy.
    # This is a general characteristic in electrochemistry.
    # The term "slower" is used in the options, which correctly describes a slow reaction.
    kinetic_conclusion = "slower"

    # --- Combine and Evaluate ---
    # Construct the correct answer based on the analysis.
    derived_answer_text = f"{thermodynamic_conclusion} – {kinetic_conclusion}"

    # The provided answer from the LLM is C.
    provided_answer_key = "C"
    options = {
        "A": "stronger – slower",
        "B": "stronger – faster",
        "C": "weaker - slower",
        "D": "weaker – faster"
    }
    llm_answer_text = options[provided_answer_key]

    # Normalize strings for a robust comparison (lowercase, remove spaces and dashes)
    def normalize(text):
        return re.sub(r'[\s–-]', '', text).lower()

    if normalize(derived_answer_text) == normalize(llm_answer_text):
        return "Correct"
    else:
        reason = []
        # Check which part is wrong
        llm_thermo, llm_kinetic = [part.strip() for part in re.split(r'[–-]', llm_answer_text)]
        
        if thermodynamic_conclusion != llm_thermo:
            reason.append(
                f"Thermodynamic part is incorrect. Oxygen is a '{thermodynamic_conclusion}' oxidant in basic solution "
                f"(since E°_basic={E_standard_basic}V < E°_acidic={E_standard_acidic}V), not '{llm_thermo}'."
            )
        
        if kinetic_conclusion != llm_kinetic:
            reason.append(
                f"Kinetic part is incorrect. Oxygen reduction is known to be kinetically '{kinetic_conclusion}', not '{llm_kinetic}'."
            )
            
        return "Incorrect. " + " ".join(reason)

# Execute the check
result = check_answer_correctness()
print(result)