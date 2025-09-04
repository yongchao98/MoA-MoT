import re

def check_answer():
    """
    Checks the correctness of the given answer for the particle physics decay problem.
    """
    # Given information from the question
    m_X = 6.0  # Mass of boson X in GeV

    # Kinematic condition: For a decay X -> f f_bar to be allowed,
    # m_X >= 2 * m_f, which means m_f <= m_X / 2.
    mass_threshold = m_X / 2.0

    # Standard Model fermion masses in GeV (approximate values)
    fermion_masses = {
        # Leptons
        'e': 0.000511,
        'mu': 0.1057,
        'tau': 1.777,
        # Quarks
        'u': 0.0022,
        'd': 0.0047,
        's': 0.095,
        'c': 1.27,
        'b': 4.18,
        't': 173.0,
    }

    # Determine the theoretically correct set of allowed decays
    theoretically_correct_decays = set()
    for fermion, mass in fermion_masses.items():
        if mass <= mass_threshold:
            theoretically_correct_decays.add(fermion)

    # The final answer provided by the LLM to be checked
    llm_answer_choice = "A"
    
    # The options as presented in the question
    options = {
        "A": "X\rightarrow c\bar{c},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}",
        "B": "X\rightarrow c\bar{c},s\bar{s},u\bar{u},d\bar{d},t\bar{t},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}",
        "C": "X\rightarrow b\bar{b},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},e^{+}e^{-}",
        "D": "X\rightarrow b\bar{b},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}"
    }

    llm_answer_string = options.get(llm_answer_choice)
    if not llm_answer_string:
        return f"Invalid answer choice '{llm_answer_choice}'. The choice must be one of {list(options.keys())}."

    # Function to parse the decay string into a set of fermion symbols
    def parse_decay_string(decay_str):
        # Remove the prefix "X\rightarrow"
        cleaned_str = decay_str.split('→')[-1]
        
        # Split by comma
        decay_parts = cleaned_str.split(',')
        
        fermion_set = set()
        # Normalize each part to a standard symbol
        for part in decay_parts:
            part = part.strip()
            if 'e' in part:
                fermion_set.add('e')
            elif 'mu' in part or 'μ' in part:
                fermion_set.add('mu')
            elif 'tau' in part or 'τ' in part:
                fermion_set.add('tau')
            elif 'c' in part:
                fermion_set.add('c')
            elif 's' in part:
                fermion_set.add('s')
            elif 'u' in part:
                fermion_set.add('u')
            elif 'd' in part:
                fermion_set.add('d')
            elif 'b' in part:
                fermion_set.add('b')
            elif 't' in part:
                fermion_set.add('t')
        return fermion_set

    # Parse the LLM's answer string
    llm_answer_decays = parse_decay_string(llm_answer_string)

    # Compare the sets
    if llm_answer_decays == theoretically_correct_decays:
        return "Correct"
    else:
        # Find the discrepancies
        missing_decays = theoretically_correct_decays - llm_answer_decays
        extra_decays = llm_answer_decays - theoretically_correct_decays
        
        error_messages = []
        if missing_decays:
            error_messages.append(f"The answer is incorrect because it omits the following kinematically allowed decays: {sorted(list(missing_decays))}.")
        
        if extra_decays:
            forbidden_details = []
            for f in extra_decays:
                forbidden_details.append(f"{f} (mass {fermion_masses[f]:.2f} GeV > {mass_threshold} GeV threshold)")
            error_messages.append(f"The answer is incorrect because it includes the following kinematically forbidden decays: {', '.join(forbidden_details)}.")
            
        return "\n".join(error_messages)

# Run the check
result = check_answer()
print(result)