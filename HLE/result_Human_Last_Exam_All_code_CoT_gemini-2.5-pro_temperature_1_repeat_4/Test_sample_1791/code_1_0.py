def solve_electron_phonon_interaction():
    """
    This function symbolically derives the effective electron-electron interaction
    mediated by phonons by integrating out the phonon fields.
    It prints the derivation step-by-step.
    """

    # Define symbolic variables as strings for clear output
    g = "g"
    q_j = "q_j"
    q_vec = "q"
    m = "m"
    w_q = "w_q"
    rho_q = "rho_q"
    rho_mq = "rho_-q"  # Represents rho with subscript -q
    H_eff = "H_eff"

    print("Derivation of the Effective Electron-Electron Interaction:\n")
    print("----------------------------------------------------------\n")

    # Step 1: Define the electron-phonon coupling vertex V_qj from the interaction term.
    # The interaction is given by g * sum_qj [ i*q_j / (2*m*w_q)^(1/2) ] * rho_q * (a_q + a_-q^dagger)
    # So, the vertex V_qj is the coefficient of rho_q * (a_q + ...).
    V_qj_expr = f"({g} * i * {q_j} / (2 * {m} * {w_q})^(1/2))"
    print(f"Step 1: The electron-phonon coupling vertex V_qj is given as:")
    print(f"V_qj = {V_qj_expr}\n")

    # Step 2: Define the vertex V_-qj. We assume w_q = w_-q.
    V_mqj_expr = f"({g} * i * (-{q_j}) / (2 * {m} * {w_q})^(1/2))"
    print(f"Step 2: The vertex for the momentum -q, V_-qj, is obtained by replacing q_j with -q_j:")
    print(f"V_-qj = {V_mqj_expr}")
    print(f"      = - V_qj\n")

    # Step 3: State the formula for the effective interaction potential V_eff(q,j) for a single mode (q,j).
    # This is obtained by completing the square in the Hamiltonian or by second-order perturbation theory.
    # The result is V_eff(q,j) = -V_qj * V_-qj / w_q
    print("Step 3: By integrating out the phonon mode (q,j), we obtain an effective interaction potential V_eff(q,j).")
    print("In the static limit, this potential is given by:")
    print(f"V_eff(q,j) = - (V_qj * V_-qj) / {w_q}")
    print(f"           = V_qj^2 / {w_q}\n")

    # Step 4: Substitute the expression for V_qj into the formula for V_eff(q,j).
    print("Step 4: Substitute the expression for V_qj into the potential formula:")
    print(f"V_eff(q,j) = {V_qj_expr}^2 / {w_q}")
    print(f"           = (-{g}^2 * {q_j}^2 / (2 * {m} * {w_q})) / {w_q}")
    V_eff_qj_final = f"(-{g}^2 * {q_j}^2 / (2 * {m} * {w_q}^2))"
    print(f"           = {V_eff_qj_final}\n")

    # Step 5: Sum over the polarization index j to get the total potential for a given q.
    print("Step 5: To get the total effective potential U_eff(q) for a given momentum q, we sum over the phonon polarizations j:")
    print(f"U_eff(q) = sum_j V_eff(q,j) = sum_j {V_eff_qj_final}")
    print(f"         = -({g}^2 / (2 * {m} * {w_q}^2)) * sum_j({q_j}^2)")
    print(f"Using the identity sum_j(q_j^2) = |{q_vec}|^2 (the squared magnitude of the vector q), we get:")
    U_eff_q_final = f"(-{g}^2 * |{q_vec}|^2 / (2 * {m} * {w_q}^2))"
    print(f"U_eff(q) = {U_eff_q_final}\n")

    # Step 6: Write down the final effective Hamiltonian term.
    print("Step 6: The final effective electron-electron interaction Hamiltonian for a given momentum q is:")
    print(f"{H_eff}(q) = U_eff(q) * {rho_q} * {rho_mq}")
    print(f"{H_eff}(q) = {U_eff_q_final} * {rho_q} * {rho_mq}\n")

    print("----------------------------------------------------------\n")
    # Final summary of the equation, showing all its components.
    print("The final equation for the effective interaction is:")
    print(f"{H_eff}(q) = - ( ( {g}^2 * |{q_vec}|^2 ) / ( 2 * {m} * {w_q}^2 ) ) * {rho_q} * {rho_mq}")

# Execute the function to print the derivation
solve_electron_phonon_interaction()