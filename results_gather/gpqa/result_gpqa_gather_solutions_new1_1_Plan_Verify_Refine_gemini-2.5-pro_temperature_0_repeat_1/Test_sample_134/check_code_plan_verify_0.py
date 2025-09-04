import re

def check_answer():
    """
    Checks the correctness of the provided answer for the particle decay question.
    """
    # --- Problem Definition ---
    # Mass of the boson X in GeV
    m_X = 6.0

    # Standard masses of fundamental fermions in GeV
    # Source: Particle Data Group (PDG) - values are consistent with the provided answers
    fermion_masses = {
        'e': 0.000511,    # electron
        'mu': 0.1057,     # muon
        'tau': 1.777,     # tau lepton
        'u': 0.0022,      # up quark
        'd': 0.0047,      # down quark
        's': 0.095,       # strange quark
        'c': 1.27,        # charm quark
        'b': 4.18,        # bottom quark
        't': 173.0,       # top quark
    }

    # The multiple-choice options provided in the question
    options_text = {
        'A': r"X\rightarrow b\bar{b},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}",
        'B': r"X\rightarrow c\bar{c},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}",
        'C': r"X\rightarrow b\bar{b},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},e^{+}e^{-}",
        'D': r"X\rightarrow c\bar{c},s\bar{s},u\bar{u},d\bar{d},t\bar{t},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}"
    }
    
    # The final answer to be checked
    final_answer_choice = 'B'

    # --- Physics Calculation ---
    # The kinematic condition for decay X -> f + f_bar is m_X >= 2 * m_f
    # This simplifies to m_f <= m_X / 2
    max_allowed_fermion_mass = m_X / 2.0

    # Determine the set of theoretically allowed decays based on the kinematic condition
    calculated_allowed_decays = set()
    for fermion, mass in fermion_masses.items():
        if mass <= max_allowed_fermion_mass:
            calculated_allowed_decays.add(fermion)

    # --- Verification ---
    # Parse the chosen answer's text to get a set of proposed decay particles
    # This handles the LaTeX formatting and extracts the base particle names
    def parse_option(text):
        particles = set()
        if r'b\bar{b}' in text: particles.add('b')
        if r's\bar{s}' in text: particles.add('s')
        if r'u\bar{u}' in text: particles.add('u')
        if r'd\bar{d}' in text: particles.add('d')
        if r'c\bar{c}' in text: particles.add('c')
        if r't\bar{t}' in text: particles.add('t')
        if r'\tau^{+}\tau^{-}' in text: particles.add('tau')
        if r'\mu^{+}\mu^{-}' in text: particles.add('mu')
        if r'e^{+}e^{-}' in text: particles.add('e')
        return particles

    proposed_decays = parse_option(options_text[final_answer_choice])

    # Compare the calculated set with the proposed set
    if calculated_allowed_decays == proposed_decays:
        return "Correct"
    else:
        # Find discrepancies to build a detailed error message
        missing_decays = calculated_allowed_decays - proposed_decays
        included_forbidden_decays = proposed_decays - calculated_allowed_decays
        
        reasons = []
        if missing_decays:
            reasons.append(f"The answer is incorrect because it is missing the following kinematically allowed decays: {sorted(list(missing_decays))}.")
        
        if included_forbidden_decays:
            reasons.append(f"The answer is incorrect because it includes the following kinematically forbidden decays: {sorted(list(included_forbidden_decays))}.")
            for f in sorted(list(included_forbidden_decays)):
                reasons.append(f"  - The decay into '{f}' is forbidden because its mass ({fermion_masses[f]} GeV) does not satisfy the condition m_f <= {max_allowed_fermion_mass} GeV.")

        return "\n".join(reasons)

# Run the check and print the result
result = check_answer()
print(result)