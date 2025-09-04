def check_particle_decay_answer():
    """
    This function checks the correctness of the provided answer for a particle decay problem.

    The problem asks for the kinematically allowed decays of a boson X with a mass of 6 GeV
    into fermion-antifermion pairs. The provided answer to check is option A.
    """

    # --- Problem Parameters ---
    # Mass of the boson X in GeV
    m_X = 6.0
    # The set of fermions listed in the given answer (Option A)
    # A) X -> c-cbar, s-sbar, u-ubar, d-dbar, tau+tau-, mu+mu-, e+e-
    answer_fermions = {'c', 's', 'u', 'd', 'tau', 'mu', 'e'}

    # --- Physics Data ---
    # A dictionary of fundamental charged fermions and their approximate masses in GeV.
    # These values are standard and sufficient for this kinematic check.
    fermion_masses = {
        # Leptons
        'e': 0.000511,   # Electron
        'mu': 0.1057,    # Muon
        'tau': 1.777,    # Tau
        # Quarks
        'u': 0.0022,     # Up
        'd': 0.0047,     # Down
        's': 0.095,      # Strange
        'c': 1.27,       # Charm
        'b': 4.18,       # Bottom
        't': 173.0       # Top
    }

    # --- Kinematic Calculation ---
    # The decay X -> f f_bar is allowed if m_X > 2 * m_f, or m_f < m_X / 2.
    kinematic_limit = m_X / 2

    # Determine the set of fermions that satisfy the kinematic constraint.
    calculated_allowed_fermions = set()
    for fermion, mass in fermion_masses.items():
        if mass < kinematic_limit:
            calculated_allowed_fermions.add(fermion)

    # --- Verification ---
    # Compare the calculated set of allowed decays with the set from the answer.
    if answer_fermions == calculated_allowed_fermions:
        return "Correct"
    else:
        # If the sets do not match, identify and report the discrepancies.
        error_messages = []

        # Check for decays incorrectly included in the answer.
        forbidden_included = answer_fermions - calculated_allowed_fermions
        for f in sorted(list(forbidden_included)):
            error_messages.append(
                f"The answer is incorrect because it includes the decay to {f}-pair. "
                f"This decay is kinematically forbidden because the fermion mass m_{f} ({fermion_masses[f]} GeV) "
                f"is not less than the kinematic limit of {kinematic_limit} GeV."
            )
        
        # Check for allowed decays that were omitted from the answer.
        allowed_omitted = calculated_allowed_fermions - answer_fermions
        for f in sorted(list(allowed_omitted)):
            error_messages.append(
                f"The answer is incorrect because it omits the allowed decay to {f}-pair. "
                f"The fermion mass m_{f} ({fermion_masses[f]} GeV) is less than the kinematic "
                f"limit of {kinematic_limit} GeV, so this decay should be included."
            )
        
        return "\n".join(error_messages)

# Execute the check and print the result.
result = check_particle_decay_answer()
print(result)