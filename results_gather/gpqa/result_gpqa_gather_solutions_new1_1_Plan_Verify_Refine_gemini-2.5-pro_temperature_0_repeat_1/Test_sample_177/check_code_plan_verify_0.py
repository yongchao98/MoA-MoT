import re

def check_qft_answer(answer_text):
    """
    Checks the correctness of the answer to the QFT problem.

    The function programmatically performs dimensional analysis to find the
    mass dimension of the coupling constant kappa and determines the
    renormalizability of the theory. It then compares this result with the
    option chosen in the provided answer.
    """

    # --- Step 1: Define the physical principles and known dimensions in 4D QFT ---
    # Mass dimension of the Lagrangian density
    dim_L = 4
    # Mass dimension of the fermion field (psi)
    # Derived from L_kin ~ psi_bar * d(psi), so 2*[psi] + 1 = 4 => [psi] = 3/2
    dim_psi = 1.5
    # Mass dimension of the field strength tensor (F_mu_nu)
    # Derived from L_kin ~ F^2, so 2*[F] = 4 => [F] = 2
    dim_F = 2.0
    # Mass dimension of sigma_mu_nu, which is made of dimensionless gamma matrices
    dim_sigma = 0.0

    # --- Step 2: Calculate the correct values based on the principles ---
    # The interaction term is L_int = kappa * psi_bar * sigma * psi * F
    # The dimensions must sum up to the dimension of the Lagrangian density.
    # [L] = [kappa] + [psi_bar] + [sigma] + [psi] + [F]
    # 4 = [kappa] + 1.5 + 0 + 1.5 + 2
    # 4 = [kappa] + 5
    correct_dim_kappa = dim_L - (dim_psi + dim_sigma + dim_psi + dim_F)

    # Determine renormalizability based on the dimension of the coupling constant.
    # A negative mass dimension implies the theory is non-renormalizable.
    correct_is_renormalizable = (correct_dim_kappa >= 0)

    # --- Step 3: Define the options from the question ---
    options = {
        'A': {'dim_kappa': 1, 'renormalizable': True, 'text': '[κ]=1, renormalizable'},
        'B': {'dim_kappa': -1, 'renormalizable': False, 'text': '[κ]=-1, not renormalizable'},
        'C': {'dim_kappa': 1, 'renormalizable': False, 'text': '[κ]=1, not renormalizable'},
        'D': {'dim_kappa': -1, 'renormalizable': True, 'text': '[κ]=-1, renormalizable'}
    }

    # --- Step 4: Parse the final answer from the provided text ---
    match = re.search(r'<<<([A-D])>>>', answer_text)
    if not match:
        return "Failure: Could not find the final answer in the format <<<X>>> in the provided text."

    chosen_option_key = match.group(1)
    chosen_option = options[chosen_option_key]

    # --- Step 5: Compare the chosen answer with the calculated correct values ---
    errors = []
    if chosen_option['dim_kappa'] != correct_dim_kappa:
        errors.append(
            f"The mass dimension of kappa is incorrect. "
            f"Expected: {correct_dim_kappa}, but the chosen answer claims: {chosen_option['dim_kappa']}."
        )

    if chosen_option['renormalizable'] != correct_is_renormalizable:
        renorm_text_expected = "renormalizable" if correct_is_renormalizable else "not renormalizable"
        renorm_text_llm = "renormalizable" if chosen_option['renormalizable'] else "not renormalizable"
        errors.append(
            f"The renormalizability is incorrect. The theory should be '{renorm_text_expected}', "
            f"but the chosen answer claims it is '{renorm_text_llm}'."
        )

    if not errors:
        return "Correct"
    else:
        # Construct a detailed error message
        correct_text = f"The correct properties are: [κ]={correct_dim_kappa}, and the theory is {'renormalizable' if correct_is_renormalizable else 'not renormalizable'}."
        chosen_text = f"The chosen option {chosen_option_key} claims: {chosen_option['text']}."
        return f"Incorrect. {' '.join(errors)} {correct_text} {chosen_text}"

# The user's provided answer text
final_answer_text = """
This problem requires a two-step analysis: first, determining the mass dimension of the coupling constant κ, and second, using that result to determine if the theory is renormalizable.

### Step 1: Dimensional Analysis of the Lagrangian Term

The fundamental principle of dimensional analysis in 4D quantum field theory (using natural units where ħ = c = 1) is that the action, S, must be dimensionless.
The action is defined as $S = \int d^4x \, \mathcal{L}$.

1.  **Dimension of the Lagrangian Density ($\mathcal{L}$):**
    *   The spacetime volume element $d^4x$ has units of (length)⁴. In natural units, mass and length have inverse dimensions, so $[M] = [L]^{-1}$.
    *   Therefore, the mass dimension of $d^4x$ is $[d^4x]_M = M^{-4}$.
    *   For the action $S$ to be dimensionless ($[S]_M = M^0$), the Lagrangian density $\mathcal{L}$ must have a mass dimension of $[\mathcal{L}]_M = M^4$.
    *   This means every individual term in the Lagrangian must have a mass dimension of 4.

2.  **Dimensions of the Fields:**
    We can find the mass dimensions of the fields by analyzing their standard kinetic terms, which must also have a dimension of $M^4$.
    *   **Fermion Field ($\psi$):** The kinetic term is $\mathcal{L}_{kin, \psi} \propto \bar{\psi}i\gamma^\mu\partial_\mu\psi$.
        *   The derivative operator $\partial_\mu$ has a mass dimension of 1 ($[\partial_\mu]_M = M^1$).
        *   The gamma matrices $\gamma^\mu$ are dimensionless.
        *   The dimensions of $\psi$ and $\bar{\psi}$ are the same.
        *   So, for the term to have dimension $M^4$, we have $[\psi]_M^2 \cdot [\partial_\mu]_M = M^4 \implies [\psi]_M^2 \cdot M^1 = M^4$.
        *   This gives $[\psi]_M^2 = M^3$, so the mass dimension of the fermion field is $[\psi]_M = 3/2$.
    *   **Field Strength Tensor ($F^{\mu\nu}$):** The kinetic term is $\mathcal{L}_{kin, A} \propto F_{\mu\nu}F^{\mu\nu}$.
        *   For this term to have dimension $M^4$, we must have $[F^{\mu\nu}]_M^2 = M^4$.
        *   This gives the mass dimension of the field strength tensor as $[F^{\mu\nu}]_M = 2$.
    *   **Sigma Tensor ($\sigma_{\mu\nu}$):** This is defined as $\sigma_{\mu\nu}=\frac{i}{2}\left[\gamma_{\mu},\gamma_{\nu}\right]$. Since the gamma matrices are dimensionless, their commutator $\sigma_{\mu\nu}$ is also dimensionless, so $[\sigma_{\mu\nu}]_M = 0$.

3.  **Dimension of the Coupling Constant ($\kappa$):**
    Now we analyze the interaction term $\mathcal{L}_{int}=\kappa\bar{\psi}\sigma_{\mu\nu}\psi F^{\mu\nu}$. Its total mass dimension must be 4.
    *   $[\mathcal{L}_{int}]_M = [\kappa]_M \cdot [\bar{\psi}]_M \cdot [\sigma_{\mu\nu}]_M \cdot [\psi]_M \cdot [F^{\mu\nu}]_M = M^4$.
    *   Substituting the dimensions we found: $[\kappa]_M \cdot M^{3/2} \cdot M^0 \cdot M^{3/2} \cdot M^2 = M^4$.
    *   $[\kappa]_M \cdot M^{(3/2 + 3/2 + 2)} = M^4$.
    *   $[\kappa]_M \cdot M^5 = M^4$.
    *   Solving for $[\kappa]_M$, we get: $[\kappa]_M = M^{-1}$.

The mass dimension of $\kappa$ is -1.

### Step 2: Renormalizability

The renormalizability of a theory can be determined by the mass dimension of its coupling constants (a method known as power-counting).
*   If a coupling constant has a **positive** mass dimension, the theory is **super-renormalizable**.
*   If a coupling constant has a **zero** mass dimension (is dimensionless), the theory is **renormalizable**.
*   If a coupling constant has a **negative** mass dimension, the theory is **non-renormalizable**.

Since the mass dimension of $\kappa$ is -1 (which is negative), the theory is **not renormalizable**.

### Conclusion

Combining our findings:
*   The mass dimension of $\kappa$ is -1.
*   The theory is not renormalizable.

This corresponds to option B. A review of the candidate answers shows that while most derive the correct physical result, many make errors in matching this result to the lettered options. The correct derivation and matching lead to B.

<<<B>>>
"""

# Run the check
result = check_qft_answer(final_answer_text)
print(result)