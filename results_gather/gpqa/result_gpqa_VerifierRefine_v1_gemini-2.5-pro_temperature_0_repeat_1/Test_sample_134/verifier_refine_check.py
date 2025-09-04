import sys
import io

def check_decay_correctness():
    """
    Checks the correctness of the provided answer for the kinematically allowed decays of boson X.
    """
    # Mass of the boson X in GeV
    m_X = 6.0

    # Kinematic condition for decay X -> f f_bar is m_X >= 2 * m_f
    # This implies m_f <= m_X / 2
    max_fermion_mass = m_X / 2.0

    # Standard Model fermion masses in GeV.
    # Using approximate values is sufficient to check the condition.
    fermion_masses = {
        # Quarks (name: mass)
        'u': 0.0022,
        'd': 0.0047,
        's': 0.095,
        'c': 1.27,
        'b': 4.18,
        't': 173.0,
        # Charged Leptons (name: mass)
        'e': 0.000511,
        'mu': 0.1057,
        'tau': 1.777,
    }

    # Determine the set of kinematically allowed decays based on the mass condition
    calculated_allowed_decays = set()
    for fermion_key, mass in fermion_masses.items():
        if mass <= max_fermion_mass:
            # Generate a standardized LaTeX string for the decay channel
            if fermion_key in ['u', 'd', 's', 'c', 'b', 't']:
                # Quark-antiquark pair, e.g., c\bar{c}
                decay_channel = f"{fermion_key}\\bar{{{fermion_key}}}"
            else:
                # Lepton-antilepton pair, e.g., e^{+}e^{-}
                lepton_symbol = {'e': 'e', 'mu': '\\mu', 'tau': '\\tau'}[fermion_key]
                decay_channel = f"{lepton_symbol}^{{+}}{lepton_symbol}^{{-}}"
            calculated_allowed_decays.add(decay_channel)

    # The decay channels listed in the proposed answer (Option C)
    # The answer is: X -> c c_bar, s s_bar, u u_bar, d d_bar, tau+ tau-, mu+ mu-, e+ e-
    # We represent this as a set of LaTeX strings for comparison.
    answer_decays = {
        "c\\bar{c}",
        "s\\bar{s}",
        "u\\bar{u}",
        "d\\bar{d}",
        "\\tau^{+}\\tau^{-}",
        "\\mu^{+}\\mu^{-}",
        "e^{+}e^{-}"
    }

    # Compare the calculated set with the answer's set
    if calculated_allowed_decays == answer_decays:
        return "Correct"
    else:
        # Find the differences to provide a specific reason for the error.
        missing_from_answer = calculated_allowed_decays - answer_decays
        extra_in_answer = answer_decays - calculated_allowed_decays

        error_messages = []
        if missing_from_answer:
            error_messages.append(f"The answer is missing the following allowed decays: {sorted(list(missing_from_answer))}.")
        
        if extra_in_answer:
            # For extra decays, explain why they are forbidden
            forbidden_reasons = []
            # Reverse map from decay string to fermion key
            fermion_map_rev = {
                'u\\bar{u}': 'u', 'd\\bar{d}': 'd', 's\\bar{s}': 's', 'c\\bar{c}': 'c', 'b\\bar{b}': 'b', 't\\bar{t}': 't',
                'e^{+}e^{-}': 'e', '\\mu^{+}\\mu^{-}': 'mu', '\\tau^{+}\\tau^{-}': 'tau'
            }
            for decay in extra_in_answer:
                fermion_key = fermion_map_rev.get(decay)
                if fermion_key:
                    mass = fermion_masses[fermion_key]
                    forbidden_reasons.append(f"{decay} (fermion mass {mass} GeV > threshold {max_fermion_mass} GeV)")
                else:
                    forbidden_reasons.append(f"{decay} (unrecognized decay)")
            error_messages.append(f"The answer incorrectly includes the following forbidden decays: {', '.join(forbidden_reasons)}.")

        return f"Incorrect. {' '.join(error_messages)}"

# Execute the check and print the result
result = check_decay_correctness()
print(result)