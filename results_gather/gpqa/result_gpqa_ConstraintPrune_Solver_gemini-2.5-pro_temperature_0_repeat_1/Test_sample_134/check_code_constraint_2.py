import sys
from io import StringIO

def check_decay_answer():
    """
    Checks the kinematically allowed decays of a boson X with a given mass.
    """
    # --- Problem Constraints ---
    mass_X = 6.0  # Mass of the boson X in GeV

    # --- Data: Fermion masses in GeV ---
    # Using standard values from the Particle Data Group (PDG)
    fermion_masses = {
        # Quarks
        'u': 0.0022,   # up quark
        'd': 0.0047,   # down quark
        's': 0.095,    # strange quark
        'c': 1.27,     # charm quark
        'b': 4.18,     # bottom quark
        't': 173.0,    # top quark
        # Charged Leptons
        'e': 0.000511, # electron
        'mu': 0.106,   # muon
        'tau': 1.777,  # tau lepton
    }

    # --- Calculation ---
    # Determine the set of kinematically allowed decays based on the mass constraint
    calculated_allowed_fermions = set()
    for fermion, mass in fermion_masses.items():
        # Kinematic condition: M_X >= 2 * m_f
        if mass_X >= 2 * mass:
            calculated_allowed_fermions.add(fermion)

    # --- Answer to Check ---
    # The answer C lists decays to c, s, u, d, tau, mu, e.
    # We represent this as a set of fermion symbols for comparison.
    # Note: 'mu' is used for muon to avoid confusion with mass 'm'.
    answer_fermions = {'c', 's', 'u', 'd', 'tau', 'mu', 'e'}

    # --- Verification ---
    # Compare the calculated set with the answer's set
    if calculated_allowed_fermions == answer_fermions:
        return "Correct"
    else:
        # If they don't match, find the discrepancies and generate an error report.
        error_messages = []

        # Check for allowed decays that are missing from the answer
        missing_decays = calculated_allowed_fermions - answer_fermions
        if missing_decays:
            error_messages.append(
                f"The answer is incorrect because it omits the following kinematically allowed decays: {sorted(list(missing_decays))}."
            )

        # Check for forbidden decays that are included in the answer
        forbidden_in_answer = answer_fermions - calculated_allowed_fermions
        if forbidden_in_answer:
            reasons = []
            for fermion in sorted(list(forbidden_in_answer)):
                mass_f = fermion_masses[fermion]
                required_mass = 2 * mass_f
                reasons.append(
                    f"X -> {fermion}-anti-{fermion} (requires M_X >= {required_mass:.3f} GeV, but M_X is only {mass_X} GeV)"
                )
            error_messages.append(
                f"The answer is incorrect because it includes the following kinematically forbidden decays: {'; '.join(reasons)}."
            )
        
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_decay_answer()
print(result)