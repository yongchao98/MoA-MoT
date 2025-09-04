import re

def check_particle_decay_answer():
    """
    Checks the correctness of the answer for the particle decay question.

    The core physics principle is that for a decay X -> f + f_bar to be
    kinematically allowed, the mass of the decaying particle (m_X) must be
    greater than or equal to the sum of the masses of the products (2 * m_f).
    Given m_X = 6 GeV, the condition is m_f <= 3 GeV.

    This code calculates the set of allowed decays based on this principle and
    compares it to the set of decays provided in the chosen answer 'C'.
    """

    # --- 1. Define problem constraints and data ---

    # Mass of the boson X in GeV
    m_X = 6.0

    # Fermion masses in GeV (using standard values from the Particle Data Group,
    # consistent with the values used in the provided answers)
    fermion_masses = {
        'u': 0.0022,    # up quark
        'd': 0.0047,    # down quark
        's': 0.095,     # strange quark
        'c': 1.27,      # charm quark
        'b': 4.18,      # bottom quark
        't': 173.0,     # top quark
        'e': 0.000511,  # electron
        'mu': 0.1057,   # muon
        'tau': 1.777,   # tau lepton
    }

    # The final answer provided in the prompt to be checked
    final_answer_key = 'C'

    # --- 2. Calculate the theoretically correct set of decays ---

    kinematic_limit = m_X / 2.0
    theoretically_allowed_decays = set()
    for fermion, mass in fermion_masses.items():
        if mass <= kinematic_limit:
            theoretically_allowed_decays.add(fermion)

    # --- 3. Parse the options from the question ---

    # Helper function to parse the decay strings into sets of particle names
    def parse_decays(decay_string):
        particles = set()
        # Split by comma, then clean up each part
        raw_particles = decay_string.split('->')[1].strip().split(',')
        for p in raw_particles:
            p_clean = p.strip()
            # Check for more specific names first to avoid ambiguity (e.g., 'u' in 'mu')
            if 'mu' in p_clean or '\\mu' in p_clean:
                particles.add('mu')
            elif 'tau' in p_clean or '\\tau' in p_clean:
                particles.add('tau')
            elif 'e' in p_clean:
                particles.add('e')
            elif 't' in p_clean:
                particles.add('t')
            elif 'b' in p_clean:
                particles.add('b')
            elif 'c' in p_clean:
                particles.add('c')
            elif 's' in p_clean:
                particles.add('s')
            elif 'd' in p_clean:
                particles.add('d')
            elif 'u' in p_clean:
                particles.add('u')
        return particles

    # The options as presented in the question
    options = {
        'A': parse_decays("X-> c-cbar,s-sbar,u-ubar,d-dbar,t-tbar,tau+tau-,mu+mu-,e+e-"),
        'B': parse_decays("X-> b-bbar,s-sbar,u-ubar,d-dbar,tau+tau-,e+e-"),
        'C': parse_decays("X-> c-cbar,s-sbar,u-ubar,d-dbar,tau+tau-,mu+mu-,e+e-"),
        'D': parse_decays("X-> b-bbar,s-sbar,u-ubar,d-dbar,tau+tau-,mu+mu-,e+e-")
    }
    
    answer_decays = options.get(final_answer_key)

    # --- 4. Compare the theoretical result with the given answer and return the verdict ---

    if answer_decays is None:
        return f"Invalid answer key '{final_answer_key}'. Key must be one of {list(options.keys())}."

    if theoretically_allowed_decays == answer_decays:
        return "Correct"
    else:
        # Find the specific errors in the answer
        forbidden_included = answer_decays - theoretically_allowed_decays
        allowed_omitted = theoretically_allowed_decays - answer_decays
        
        reasons = []
        if forbidden_included:
            # Check the mass of the forbidden particles to provide a clear reason
            for particle in sorted(list(forbidden_included)):
                reasons.append(f"the answer incorrectly includes the forbidden decay to '{particle}' (mass {fermion_masses[particle]} GeV > {kinematic_limit} GeV)")
        if allowed_omitted:
            reasons.append(f"the answer omits the allowed decay(s) to {sorted(list(allowed_omitted))}")
            
        return f"Incorrect. The kinematic condition is not satisfied. Specifically, {'; '.join(reasons)}."

# Execute the check and print the result
print(check_particle_decay_answer())