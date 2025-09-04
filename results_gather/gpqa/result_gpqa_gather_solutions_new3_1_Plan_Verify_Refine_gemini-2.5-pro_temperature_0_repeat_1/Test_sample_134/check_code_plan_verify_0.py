import re

def check_correctness_of_particle_decay():
    """
    Checks the correctness of the proposed answer for the kinematically allowed decays of boson X.

    The function performs the following steps:
    1. Defines the physical parameters: mass of boson X and masses of all fundamental fermions.
    2. Calculates the theoretically correct set of allowed decays based on the kinematic condition: m_f <= m_X / 2.
    3. Parses the proposed answer ('<<<B>>>') to identify the chosen option.
    4. Maps the chosen option to its corresponding set of decay products.
    5. Compares the theoretically correct set with the set from the proposed answer.
    6. Returns "Correct" if they match, or a detailed error message explaining the discrepancy.
    """
    # --- Step 1: Define physical constraints and data ---
    m_X = 6.0  # Mass of boson X in GeV
    mass_threshold = m_X / 2  # Maximum allowed mass for a fermion decay product

    # Fermion masses in GeV (using values from the provided analysis for consistency)
    fermion_masses = {
        'e': 0.000511,   # electron
        'mu': 0.1057,    # muon
        'tau': 1.777,    # tau lepton
        'u': 0.0022,     # up quark
        'd': 0.0047,     # down quark
        's': 0.095,      # strange quark
        'c': 1.27,       # charm quark
        'b': 4.18,       # bottom quark
        't': 173.0,      # top quark
    }

    # --- Step 2: Determine the theoretically correct set of allowed decays ---
    correct_allowed_fermions = set()
    for fermion, mass in fermion_masses.items():
        if mass <= mass_threshold:
            correct_allowed_fermions.add(fermion)

    # --- Step 3: Parse the proposed answer ---
    # The final answer provided by the LLM analysis to be checked
    proposed_answer_text = "<<<B>>>"

    # Define the decay channels for each option as listed in the original question
    options = {
        'A': {'c', 's', 'u', 'd', 't', 'tau', 'mu', 'e'},
        'B': {'c', 's', 'u', 'd', 'tau', 'mu', 'e'},
        'C': {'b', 's', 'u', 'd', 'tau', 'e'},
        'D': {'b', 's', 'u', 'd', 'tau', 'mu', 'e'}
    }

    match = re.search(r'<<<([A-D])>>>', proposed_answer_text)
    if not match:
        return "Error: Invalid answer format. The answer should be in the format <<<X>>> where X is A, B, C, or D."

    proposed_option_key = match.group(1)
    proposed_fermions = options.get(proposed_option_key)

    # --- Step 4: Compare the correct set with the proposed set ---
    if correct_allowed_fermions == proposed_fermions:
        return "Correct"
    else:
        # Analyze the discrepancies to provide a clear reason for the error
        missing_decays = correct_allowed_fermions - proposed_fermions
        extra_decays = proposed_fermions - correct_allowed_fermions

        error_messages = []
        if extra_decays:
            # These are decays that are kinematically forbidden but were included in the answer
            reasons = [f"{f} (mass {fermion_masses[f]} GeV > {mass_threshold} GeV)" for f in sorted(list(extra_decays))]
            error_messages.append(f"The answer incorrectly includes forbidden decays: {', '.join(reasons)}.")

        if missing_decays:
            # These are decays that are kinematically allowed but were omitted from the answer
            reasons = [f"{f} (mass {fermion_masses[f]} GeV <= {mass_threshold} GeV)" for f in sorted(list(missing_decays))]
            error_messages.append(f"The answer omits allowed decays: {', '.join(reasons)}.")

        return "Incorrect. " + " ".join(error_messages)

# Run the check and print the result
result = check_correctness_of_particle_decay()
print(result)