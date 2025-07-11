def derive_effective_interaction():
    """
    This function presents a symbolic derivation of the effective
    electron-electron interaction mediated by phonons.
    """
    print("Derivation of the Phonon-Mediated Effective Electron-Electron Interaction")
    print("="*70)
    print("Step 1: Define the Hamiltonians")
    print("-" * 70)
    H_ph = "H_ph = sum_{q, j} w_q * a_{q,j}^dagger * a_{q,j}"
    H_el_ph_full = "g * sum_{k, q, j} (i * q_j / (2 * m * w_q)^(1/2)) * n_{q} * (a_{q,j} + a_{-q,j}^dagger)"
    H_el_ph = "H_el-ph = sum_{q, j} M_{q,j} * n_q * (a_{q,j} + a_{-q,j}^dagger)"
    M_qj = "g * i * q_j / (2 * m * w_q)^(1/2)"
    print(f"The given phonon Hamiltonian is: {H_ph}")
    print(f"The given electron-phonon interaction is: {H_el_ph_full}")
    print("We can write the interaction part more compactly as:")
    print(f"H_el-ph = {H_el_ph}")
    print(f"where the coupling vertex M_{q,j} is: {M_qj}")
    print("Here, n_q is the electron density operator, and a, a^dagger are phonon operators.")
    print("\nOur goal is to integrate out the phonon fields (a, a^dagger) to find an effective interaction between electron densities (n_q).")
    print()

    print("Step 2: The Effective Interaction via Path Integral")
    print("-" * 70)
    print("Integrating out the quadratic phonon fields in the path integral is equivalent to")
    print("calculating the second-order perturbation contribution from H_el-ph.")
    print("The resulting effective interaction potential U_eff(q, i*nu_n) between electron densities")
    print("n_q and n_{-q} is given by the following expression in Matsubara frequency space:")
    U_eff_formula = "-sum_j M_{q,j} * M_{-q,j} * D_{q,j}(i*nu_n)"
    print(f"\nU_eff(q, i*nu_n) = {U_eff_formula}\n")
    print("where D_{q,j}(i*nu_n) is the bare phonon propagator.")
    print()

    print("Step 3: Evaluate the Components")
    print("-" * 70)
    
    # Component 1: The product of coupling vertices
    print("First, we compute the product of the vertices M_{q,j} and M_{-q,j}:")
    print(f"M_{q,j} = {M_qj}")
    M_minus_qj = "g * i * (-q_j) / (2 * m * w_q)^(1/2)"
    print(f"M_{-q,j} = {M_minus_qj} (since w_{-q} = w_q)")
    Product_M = "(g * i * q_j) * (g * i * -q_j) / (2 * m * w_q)"
    Product_M_simplified = "g^2 * q_j^2 / (2 * m * w_q)"
    print(f"M_{q,j} * M_{-q,j} = {Product_M} = {Product_M_simplified}")
    print()

    # Component 2: The phonon propagator
    print("Second, we need the phonon propagator D_{q,j}(i*nu_n). It is the Fourier transform of:")
    Propagator_def = "<T_tau (a_{q,j}(tau) + a_{-q,j}^dagger(tau)) * (a_{-q,j}(0) + a_{q,j}^dagger(0))>"
    print(f"D_{q,j}(tau) = {Propagator_def}")
    print("Its Fourier transform to Matsubara frequencies (nu_n) is well-known:")
    Propagator_fourier = "2 * w_q / (nu_n^2 + w_q^2)"
    print(f"D_{q,j}(i*nu_n) = {Propagator_fourier}")
    print()

    print("Step 4: Assemble the Final Expression")
    print("-" * 70)
    print("Now we substitute these parts back into the formula for U_eff(q, i*nu_n):")
    U_eff_intermediate = f"-sum_j ({Product_M_simplified}) * ({Propagator_fourier})"
    print(f"U_eff(q, i*nu_n) = {U_eff_intermediate}")
    
    print("\nSimplifying the expression inside the sum:")
    term_j = "-g^2 * q_j^2 / (m * (nu_n^2 + w_q^2))"
    print(f"The term for a given polarization j is: {term_j}")
    
    print("\nFinally, we perform the sum over the polarization index j (e.g., j=x,y,z in 3D).")
    sum_qj2 = "sum_j q_j^2 = q^2"
    print(f"We use the identity: {sum_qj2}, where q is the magnitude of the momentum vector q.")
    
    print("\nThis yields the final result for the effective electron-electron interaction potential.")
    print("\n--- FINAL RESULT ---")
    print("The effective interaction potential U_eff for a momentum transfer q is:")
    
    # Printing the structure of the final equation piece by piece to show "each number".
    final_expression = "U_eff(q, i*nu_n) = -(g^2 * q^2) / (m * (nu_n^2 + w_q^2))"
    print(final_expression)
    
    print("\nBreaking down the final equation into its components:")
    numerator_constant = -1
    g_exponent = 2
    q_exponent = 2
    print(f"Numerator = ({numerator_constant}) * g^{g_exponent} * q^{q_exponent}")

    m_exponent = 1
    nu_n_exponent = 2
    w_q_exponent = 2
    print(f"Denominator = m^{m_exponent} * (nu_n^{nu_n_exponent} + w_q^{w_q_exponent})")


if __name__ == '__main__':
    derive_effective_interaction()