import sys
import io

def check_decay_correctness():
    """
    Checks the kinematically allowed decays for a boson X with a mass of 6 GeV.

    The kinematic condition for a decay X -> f f_bar is m_X > 2 * m_f.
    Given m_X = 6 GeV, the condition becomes m_f < 3 GeV.

    This function defines the masses of fundamental fermions, calculates the set of
    allowed decay products, and compares it against the set provided in the
    selected answer (Option C).
    """
    # Mass of the boson X in GeV
    m_X = 6.0

    # Masses of fundamental fermions in GeV
    # Source: Particle Data Group (PDG) values, consistent with the LLM's reasoning
    fermion_masses = {
        # Leptons
        'e': 0.000511,    # Electron
        'mu': 0.1057,     # Muon
        'tau': 1.777,     # Tau
        # Quarks
        'u': 0.0022,      # Up
        'd': 0.0047,      # Down
        's': 0.095,       # Strange
        'c': 1.27,        # Charm
        'b': 4.18,        # Bottom
        't': 173.0,       # Top
    }

    # The options provided in the question, represented as sets of fermion symbols.
    # Note: The decay is to a fermion-antifermion pair, so we just list the fermion.
    options = {
        'A': {'b', 's', 'u', 'd', 'tau', 'mu', 'e'},
        'B': {'c', 's', 'u', 'd', 't', 'tau', 'mu', 'e'},
        'C': {'c', 's', 'u', 'd', 'tau', 'mu', 'e'},
        'D': {'b', 's', 'u', 'd', 'tau', 'e'},
    }

    # The answer provided by the LLM
    llm_answer_key = 'C'

    # Kinematic condition: The fermion mass must be less than half the boson's mass.
    max_allowed_fermion_mass = m_X / 2.0

    # Determine the set of kinematically allowed fermions based on the condition
    calculated_allowed_fermions = set()
    for fermion, mass in fermion_masses.items():
        if mass < max_allowed_fermion_mass:
            calculated_allowed_fermions.add(fermion)

    # Get the set of fermions from the LLM's chosen answer
    llm_answer_fermions = options.get(llm_answer_key)

    if llm_answer_fermions is None:
        return f"The provided answer key '{llm_answer_key}' is not a valid option (A, B, C, or D)."

    # Compare the calculated set with the answer's set
    if calculated_allowed_fermions == llm_answer_fermions:
        return "Correct"
    else:
        # Find the discrepancies to provide a detailed reason for the error
        missing_from_answer = calculated_allowed_fermions - llm_answer_fermions
        extra_in_answer = llm_answer_fermions - calculated_allowed_fermions

        error_messages = []
        if extra_in_answer:
            details = []
            for f in sorted(list(extra_in_answer)):
                details.append(f"'{f}' (mass ~{fermion_masses[f]:.2f} GeV)")
            error_messages.append(
                f"The answer is incorrect because it includes forbidden decays. "
                f"The following particles are too heavy to be produced: {', '.join(details)}. "
                f"The kinematic condition is m_fermion < {max_allowed_fermion_mass} GeV."
            )
        
        if missing_from_answer:
            details = []
            for f in sorted(list(missing_from_answer)):
                details.append(f"'{f}' (mass ~{fermion_masses[f]:.2f} GeV)")
            error_messages.append(
                f"The answer is incorrect because it omits allowed decays. "
                f"The following particles should be included: {', '.join(details)}. "
                f"Their masses are less than the {max_allowed_fermion_mass} GeV threshold."
            )
            
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_decay_correctness()
print(result)