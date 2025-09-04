import re

def check_particle_decay_answer():
    """
    Checks the correctness of the LLM's answer regarding kinematically allowed decays.

    The function performs the following steps:
    1. Defines the physical constants (fermion masses) and problem parameters (boson mass).
    2. Calculates the correct set of allowed decays based on the kinematic condition: m_boson >= 2 * m_fermion.
    3. Parses the multiple-choice options as presented in the original question.
    4. Extracts the final choice (e.g., 'B') from the provided LLM answer text.
    5. Compares the set of decays corresponding to the LLM's choice with the correctly calculated set.
    6. Returns "Correct" if they match, or a detailed error message if they do not.
    """
    # 1. Define physical constants and problem parameters
    fermion_masses = {
        # Leptons (in GeV/c^2)
        'e': 0.000511,
        'mu': 0.1057,
        'tau': 1.777,
        # Quarks (in GeV/c^2)
        'u': 0.0022,
        'd': 0.0047,
        's': 0.095,
        'c': 1.27,
        'b': 4.18,
        't': 173.0,
    }
    boson_mass = 6.0  # in GeV

    # 2. Determine the correct set of allowed decays based on physics
    mass_threshold = boson_mass / 2
    correct_allowed_fermions = {
        fermion for fermion, mass in fermion_masses.items() if mass <= mass_threshold
    }

    # 3. Define the options as presented in the question prompt
    # A) X->b,s,u,d,tau,mu,e
    # B) X->c,s,u,d,tau,mu,e
    # C) X->c,s,u,d,t,tau,mu,e
    # D) X->b,s,u,d,tau,e
    options = {
        'A': {'b', 's', 'u', 'd', 'tau', 'mu', 'e'},
        'B': {'c', 's', 'u', 'd', 'tau', 'mu', 'e'},
        'C': {'c', 's', 'u', 'd', 't', 'tau', 'mu', 'e'},
        'D': {'b', 's', 'u', 'd', 'tau', 'e'}
    }

    # 4. Parse the LLM's final answer from the provided text
    llm_answer_text = """
    Here is a step-by-step analysis of the problem to determine the correct final answer.

    ### Step-by-step Analysis

    1.  **Understand the Physical Principle:** The problem asks for the "kinematically allowed" decays of a boson X. In particle physics, a decay is kinematically allowed if it conserves energy. For a particle at rest decaying into other particles, this simplifies to the condition that the mass of the initial particle must be greater than or equal to the sum of the masses of the final particles.

    2.  **Formulate the Kinematic Condition:** The interaction Lagrangian, $\mathcal{L}\left(y\right)=-\lambda{f}\bar{\psi}{f}\left(y\right)\left(\upsilon+X\left(y\right)\right)\psi_{f}\left(y\right)$, contains the term $-\lambda_f X \bar{\psi}_f \psi_f$. This describes the decay of the boson X into a fermion-antifermion pair ($f\bar{f}$). The kinematic condition for this decay, $X \rightarrow f\bar{f}$, is:
        $$ m_X \ge m_f + m_{\bar{f}} $$
        Since a fermion ($f$) and its corresponding antifermion ($\bar{f}$) have identical masses, the condition becomes:
        $$ m_X \ge 2m_f $$

    3.  **Apply the Given Values:** We are given the mass of the boson X as $m_X = 6$ GeV. Plugging this into the condition gives the requirement for the fermion's mass, $m_f$:
        $$ 6 \text{ GeV} \ge 2m_f $$
        Solving for $m_f$, we find the condition for an allowed decay:
        $$ m_f \le 3 \text{ GeV} $$
        This means the boson X can decay into any fermion-antifermion pair as long as the individual fermion's mass is less than or equal to 3 GeV.

    4.  **Check Fermion Masses:** We must now check this condition for all fundamental fermions using their standard approximate masses.

        *   **Leptons:**
            *   Electron ($e$): $m_e \approx 0.000511$ GeV. Since $0.000511 \le 3$, the decay $X \rightarrow e^{+}e^{-}$ is **allowed**.
            *   Muon ($\mu$): $m_\mu \approx 0.1057$ GeV. Since $0.1057 \le 3$, the decay $X \rightarrow \mu^{+}\mu^{-}$ is **allowed**.
            *   Tau ($\tau$): $m_\tau \approx 1.777$ GeV. Since $1.777 \le 3$, the decay $X \rightarrow \tau^{+}\tau^{-}$ is **allowed**.

        *   **Quarks:**
            *   Up ($u$): $m_u \approx 0.0022$ GeV. Since $0.0022 \le 3$, the decay $X \rightarrow u\bar{u}$ is **allowed**.
            *   Down ($d$): $m_d \approx 0.0047$ GeV. Since $0.0047 \le 3$, the decay $X \rightarrow d\bar{d}$ is **allowed**.
            *   Strange ($s$): $m_s \approx 0.095$ GeV. Since $0.095 \le 3$, the decay $X \rightarrow s\bar{s}$ is **allowed**.
            *   Charm ($c$): $m_c \approx 1.27$ GeV. Since $1.27 \le 3$, the decay $X \rightarrow c\bar{c}$ is **allowed**.
            *   Bottom ($b$): $m_b \approx 4.18$ GeV. Since $4.18 > 3$, the decay $X \rightarrow b\bar{b}$ is **forbidden**.
            *   Top ($t$): $m_t \approx 173$ GeV. Since $173 > 3$, the decay $X \rightarrow t\bar{t}$ is **forbidden**.

    5.  **Compile the Final List and Compare with Options:**
        The complete set of kinematically allowed decays is:
        $$ X\rightarrow c\bar{c},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-} $$
        Now we must carefully compare this correct list with the options provided in the question prompt:
        *   A) $X\rightarrow b\bar{b},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}$: **Incorrect**. Includes the forbidden $b\bar{b}$ decay and omits the allowed $c\bar{c}$ decay.
        *   B) $X\rightarrow c\bar{c},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}$: **Correct**. This list perfectly matches all the allowed decays and excludes all forbidden ones.
        *   C) $X\rightarrow c\bar{c},s\bar{s},u\bar{u},d\bar{d},t\bar{t},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}$: **Incorrect**. Includes the forbidden $t\bar{t}$ decay.
        *   D) $X\rightarrow b\bar{b},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},e^{+}e^{-}$: **Incorrect**. Includes the forbidden $b\bar{b}$ decay and omits the allowed $c\bar{c}$ and $\mu^{+}\mu^{-}$ decays.

    ### Analysis of Candidate Answers and Careful Points

    *   **Consensus on Physics:** The provided candidate answers are nearly unanimous in their physical reasoning. They all correctly identify the kinematic condition ($m_X \ge 2m_f$), calculate the fermion mass threshold ($m_f \le 3$ GeV), and derive the same correct set of allowed decay products: $\{c, s, u, d, \tau, \mu, e\}$.
    *   **Discrepancy in Final Choice:** The major issue is a widespread inconsistency in matching this correct set of decays to the lettered options (A, B, C, D). For example, many candidates who chose "A" or "C" correctly described the contents of option "B" in their reasoning. This strongly suggests that the agents were working with different or scrambled versions of the multiple-choice options.
    *   **Conclusion:** The physical analysis robustly points to a single set of allowed decays. The errors in the candidate answers are clerical (matching the result to the wrong letter), not physical. Based on a direct and careful comparison of the derived result with the options as stated in the problem prompt, option B is the only correct choice.

    <<<B>>>
    """
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a final answer in the format <<<A>>>."
    llm_choice = match.group(1)

    # 5. Compare the LLM's choice with the correct answer
    llm_selected_fermions = options.get(llm_choice)

    if llm_selected_fermions == correct_allowed_fermions:
        return "Correct"
    else:
        # Analyze the error for a detailed message
        incorrectly_included = llm_selected_fermions - correct_allowed_fermions
        incorrectly_omitted = correct_allowed_fermions - llm_selected_fermions
        
        error_messages = []
        if incorrectly_included:
            details = []
            for f in sorted(list(incorrectly_included)):
                details.append(f"'{f}' (mass {fermion_masses[f]} GeV > {mass_threshold} GeV)")
            error_messages.append(f"The answer incorrectly includes forbidden decays to: {', '.join(details)}.")
        
        if incorrectly_omitted:
            details = []
            for f in sorted(list(incorrectly_omitted)):
                details.append(f"'{f}' (mass {fermion_masses[f]} GeV <= {mass_threshold} GeV)")
            error_messages.append(f"The answer incorrectly omits allowed decays to: {', '.join(details)}.")
            
        return (f"Incorrect. The chosen option '{llm_choice}' is wrong. "
                f"The correct set of allowed decays is {sorted(list(correct_allowed_fermions))}. "
                + " ".join(error_messages))

# Execute the check and print the result
result = check_particle_decay_answer()
print(result)