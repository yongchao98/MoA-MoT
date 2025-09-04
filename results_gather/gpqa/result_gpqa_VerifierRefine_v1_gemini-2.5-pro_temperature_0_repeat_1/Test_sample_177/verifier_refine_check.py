import sys
from io import StringIO

def check_qft_answer():
    """
    This function programmatically checks the dimensional analysis of a QFT interaction term.
    It calculates the mass dimension of the coupling constant and determines the theory's
    renormalizability, then compares it to the provided answer.
    """
    
    # --- Step 1: Define fundamental principles in 4D spacetime ---
    # In natural units (hbar=c=1), mass dimension [M] is used.
    # [Length] = [Time] = -1.
    # The action S = integral(d^4x L) must be dimensionless, [S] = 0.
    spacetime_dimension = 4
    
    # The mass dimension of the spacetime volume element d^4x is [d^4x] = -4.
    dim_d4x = -spacetime_dimension
    
    # For S to be dimensionless, the Lagrangian density L must have [L] = 4.
    dim_L = -dim_d4x
    
    # --- Step 2: Determine mass dimensions of the fields from their kinetic terms ---
    
    # Fermion field psi: The kinetic term is L_kin ~ psi_bar * (i * gamma^mu * d_mu) * psi.
    # This term must have dimension [L] = 4.
    # [psi_bar] + [d_mu] + [psi] = 4.
    # Gamma matrices are dimensionless. [d_mu] = 1. [psi_bar] = [psi].
    # 2 * [psi] + 1 = 4  =>  [psi] = 3/2.
    dim_partial_derivative = 1
    dim_psi = (dim_L - dim_partial_derivative) / 2.0
    
    # Gauge field strength tensor F^munu: The kinetic term is L_kin ~ -1/4 * F_munu * F^munu.
    # This term must have dimension [L] = 4.
    # 2 * [F^munu] = 4  =>  [F^munu] = 2.
    dim_F = dim_L / 2.0
    
    # Sigma tensor sigma_munu = i/2 * [gamma_mu, gamma_nu].
    # Since gamma matrices are dimensionless constant matrices, their commutator is also dimensionless.
    dim_sigma = 0
    
    # --- Step 3: Calculate the mass dimension of the coupling constant kappa ---
    
    # The interaction term L_int = kappa * psi_bar * sigma_munu * psi * F^munu must also have dimension [L] = 4.
    # [kappa] + [psi_bar] + [sigma_munu] + [psi] + [F^munu] = 4.
    # Note: [psi_bar] = [psi].
    # [kappa] + 2 * [psi] + [sigma_munu] + [F^munu] = 4.
    
    # Substitute the calculated dimensions:
    # [kappa] + 2 * (3/2) + 0 + 2 = 4
    # [kappa] + 3 + 2 = 4
    # [kappa] + 5 = 4
    # [kappa] = -1
    
    calculated_dim_kappa = dim_L - (2 * dim_psi + dim_sigma + dim_F)
    
    # --- Step 4: Determine if the theory is renormalizable ---
    
    # The renormalizability is determined by the mass dimension of the coupling constant.
    # [g] > 0  => Super-renormalizable
    # [g] = 0  => Renormalizable
    # [g] < 0  => Non-renormalizable
    
    if calculated_dim_kappa < 0:
        calculated_renormalizability = "not renormalizable"
    elif calculated_dim_kappa == 0:
        calculated_renormalizability = "renormalizable"
    else:
        calculated_renormalizability = "super-renormalizable"
        
    # --- Step 5: Compare with the given answer (Option C) ---
    # Option C states: The mass dimension [kappa] = -1. The theory is not renormalizable.
    
    expected_dim_kappa = -1.0
    expected_renormalizability = "not renormalizable"
    
    errors = []
    
    # Check the mass dimension
    if calculated_dim_kappa != expected_dim_kappa:
        errors.append(
            f"Dimension Mismatch: The calculated mass dimension of kappa is [{calculated_dim_kappa}], "
            f"but the answer expects [{expected_dim_kappa}]."
        )
        
    # Check the renormalizability
    if calculated_renormalizability != expected_renormalizability:
        errors.append(
            f"Renormalizability Mismatch: Based on the calculated dimension [{calculated_dim_kappa}], "
            f"the theory is '{calculated_renormalizability}', but the answer states it is "
            f"'{expected_renormalizability}'."
        )
        
    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        # Combine all found errors into a single message.
        return "Incorrect. " + " ".join(errors)

# Execute the check and print the result.
result = check_qft_answer()
print(result)