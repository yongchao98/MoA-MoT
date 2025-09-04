import re

def check_electrochemistry_answer():
    """
    This function checks the correctness of the provided answer to the electrochemistry question.
    It verifies the thermodynamic and kinetic properties of oxygen based on established principles.
    """
    
    # --- Part 1: Define the ground truth based on electrochemistry principles ---

    # 1a. Thermodynamic Analysis
    # The strength of an oxidant is determined by its standard reduction potential (E°).
    # A higher (more positive) E° indicates a stronger oxidant.
    E_potential_acidic = 1.23  # V for O₂ + 4H⁺ + 4e⁻ → 2H₂O
    E_potential_basic = 0.40   # V for O₂ + 2H₂O + 4e⁻ → 4OH⁻

    # The question asks about oxygen in basic solutions. We compare its potential to the acidic case.
    if E_potential_basic < E_potential_acidic:
        correct_thermo_term = "weaker"
    else:
        correct_thermo_term = "stronger"

    # 1b. Kinetic Analysis
    # The oxygen reduction reaction (ORR) is famously kinetically "sluggish" or slow,
    # especially in acidic media, due to a high activation energy.
    # Compared to alkaline (basic) media, the reaction in acidic media is generally slower.
    correct_kinetic_term = "slower"

    # 1c. Form the correct answer phrase
    # The question format is "weaker/stronger – faster/slower"
    correct_phrase = f"{correct_thermo_term} - {correct_kinetic_term}"

    # --- Part 2: Define the problem's options and the given answer ---

    # The options provided in the question
    options = {
        "A": "stronger – slower",
        "B": "weaker – faster",
        "C": "weaker - slower",
        "D": "stronger – faster"
    }
    
    # Standardize dashes in options for robust comparison
    standardized_options = {k: v.replace(" – ", " - ").replace("–", "-") for k, v in options.items()}

    # The answer provided by the LLM to be checked
    llm_answer_text = "<<<C>>>"
    
    # --- Part 3: Verification Logic ---

    # Find which option letter corresponds to the correct phrase
    correct_option_letter = None
    for letter, phrase in standardized_options.items():
        if phrase == correct_phrase:
            correct_option_letter = letter
            break
            
    if correct_option_letter is None:
        return "Checker Error: The derived correct phrase '{correct_phrase}' does not match any of the provided options."

    # Extract the letter from the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format. The provided answer was '{llm_answer_text}', but it should be in the format '<<<X>>>' where X is A, B, C, or D."
        
    llm_option_letter = match.group(1)

    # Compare the LLM's answer with the derived correct answer
    if llm_option_letter == correct_option_letter:
        return "Correct"
    else:
        reason = (
            f"The provided answer is '{llm_option_letter}', but the correct answer is '{correct_option_letter}'.\n"
            f"Reasoning:\n"
            f"1. Thermodynamics: The standard reduction potential of oxygen in basic solution (E° = {E_potential_basic}V) is lower than in acidic solution (E° = {E_potential_acidic}V). A lower potential means it is a 'weaker' oxidant. The first term must be 'weaker'.\n"
            f"2. Kinetics: The oxygen reduction reaction is known to be kinetically 'slower' in acidic solutions due to a high activation energy and sluggish multi-electron transfer process.\n"
            f"Therefore, the correct combination is '{correct_phrase}', which corresponds to option {correct_option_letter}."
        )
        return reason

# Execute the check and print the result
result = check_electrochemistry_answer()
print(result)