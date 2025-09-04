import sys
import io

def check_qft_answer(llm_answer_choice: str):
    """
    Checks the correctness of the answer to the QFT problem.

    The function calculates the mass dimension of the coupling constant kappa
    and determines the renormalizability of the theory based on first principles.
    It then compares this result with the provided answer choice.

    Args:
        llm_answer_choice (str): The letter corresponding to the chosen answer (e.g., 'A', 'B').

    Returns:
        str: "Correct" if the answer is right, otherwise a string explaining the error.
    """
    # --- Step 1: Define the problem's constraints and knowns ---
    # In natural units (hbar=c=1) and 4D spacetime:
    # The action S = integral(L d^4x) is dimensionless.
    # Mass dimension of spacetime volume [d^4x] = -4.
    # Therefore, mass dimension of Lagrangian density [L] must be 4.
    
    # Standard mass dimensions of fields from their kinetic terms:
    # L_fermion = psi_bar * (i * gamma^mu * d_mu - m) * psi => [psi] = 3/2
    # L_gauge = -1/4 * F_munu * F^munu => [F_munu] = 2
    
    dims = {
        'L_int': 4,
        'psi': 3/2,
        'F_munu': 2,
        'sigma_munu': 0,  # sigma_munu is built from dimensionless gamma matrices
    }

    # The interaction Lagrangian is L_int = kappa * psi_bar * sigma_munu * psi * F^munu
    # The sum of dimensions of terms in L_int must equal dim['L_int']
    # [kappa] + [psi_bar] + [sigma_munu] + [psi] + [F_munu] = 4
    # Note: [psi_bar] = [psi]

    # --- Step 2: Calculate the mass dimension of kappa ---
    calculated_kappa_dim = dims['L_int'] - (dims['psi'] + dims['sigma_munu'] + dims['psi'] + dims['F_munu'])

    # --- Step 3: Determine renormalizability based on the calculated dimension ---
    # A theory is non-renormalizable by power-counting if any coupling constant
    # has a negative mass dimension.
    is_theory_renormalizable = (calculated_kappa_dim >= 0)

    # --- Step 4: Define the options and parse the given answer ---
    options = {
        'A': {'kappa_dim': -1, 'renormalizable': True, 'text': 'The mass dimension [kappa]=-1. The theory is renormalizable.'},
        'B': {'kappa_dim': -1, 'renormalizable': False, 'text': 'The mass dimension [kappa]=-1. The theory is not renormalizable.'},
        'C': {'kappa_dim': 1, 'renormalizable': True, 'text': 'The mass dimension [kappa]=1. The theory is renormalizable.'},
        'D': {'kappa_dim': 1, 'renormalizable': False, 'text': 'The mass dimension [kappa]=1. The theory is not renormalizable.'}
    }

    if llm_answer_choice.upper() not in options:
        return f"Error: The provided answer '{llm_answer_choice}' is not a valid option (A, B, C, or D)."

    chosen_option = options[llm_answer_choice.upper()]
    
    # --- Step 5: Compare the calculated results with the chosen answer ---
    errors = []
    
    # Check 1: Mass dimension of kappa
    if calculated_kappa_dim != chosen_option['kappa_dim']:
        error_msg = (
            f"Constraint Violated: Mass dimension of the coupling constant kappa.\n"
            f"Reason: The answer claims [kappa] = {chosen_option['kappa_dim']}, but the correct value is {calculated_kappa_dim}.\n"
            f"Derivation: The Lagrangian density L_int must have a mass dimension of 4. "
            f"From L_int = k * psi_bar * sigma * psi * F, we get [k] + [psi_bar] + [sigma] + [psi] + [F] = 4. "
            f"With [psi]=3/2, [sigma]=0, and [F]=2, the equation is [k] + 3/2 + 0 + 3/2 + 2 = 4, which solves to [k] = -1."
        )
        errors.append(error_msg)

    # Check 2: Renormalizability
    if is_theory_renormalizable != chosen_option['renormalizable']:
        renorm_text_correct = "renormalizable" if is_theory_renormalizable else "not renormalizable"
        renorm_text_claimed = "renormalizable" if chosen_option['renormalizable'] else "not renormalizable"
        error_msg = (
            f"Constraint Violated: Renormalizability of the theory.\n"
            f"Reason: The answer claims the theory is {renorm_text_claimed}, but it is {renorm_text_correct}.\n"
            f"Derivation: A theory is considered non-renormalizable by power-counting if its coupling constant has a negative mass dimension. "
            f"Since the calculated [kappa] = {calculated_kappa_dim} is negative, the theory is not renormalizable."
        )
        errors.append(error_msg)

    # --- Step 6: Return the final verdict ---
    if not errors:
        return "Correct"
    else:
        return "\n\n".join(errors)

# The provided answer from the other LLM is 'B'.
# Let's run the check.
llm_answer = 'B' 
result = check_qft_answer(llm_answer)
print(result)