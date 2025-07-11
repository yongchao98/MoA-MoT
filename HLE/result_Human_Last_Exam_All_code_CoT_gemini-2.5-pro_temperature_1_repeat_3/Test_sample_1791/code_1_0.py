def display_effective_interaction():
    """
    This function prints the derived effective electron-electron interaction potential.
    Each variable in the equation is explicitly named and printed.
    """
    
    # Define the components of the equation as strings
    potential = "V_eff(q, i*v_n)"
    sum_term = "sum_j"
    g_sq = "g^2"
    q_j_sq = "q_j^2"
    m = "m"
    v_n_sq = "v_n^2"
    w_q_sq = "w_q^2"
    
    # Assemble the final equation string
    # V_eff(q, i*v_n) = sum_j (g^2 * q_j^2 / m) * (1 / (v_n^2 + w_q^2))
    equation = f"{potential} = {sum_term} ({g_sq} * {q_j_sq} / {m}) * (1 / ({v_n_sq} + {w_q_sq}))"

    print("The effective electron-electron interaction potential V_eff for a given momentum transfer q is:")
    print("-" * 60)
    
    # Print the final equation with each part identified
    print(f"The potential V_eff as a function of momentum q and Matsubara frequency v_n is:")
    print(f"  {potential} = {sum_term} (({g_sq})*({q_j_sq}) / {m}) * (1 / (({v_n_sq}) + ({w_q_sq})))")
    
    print("\nWhere:")
    print(f"  V_eff: The effective interaction potential.")
    print(f"  q: The momentum transfer vector.")
    print(f"  j: The phonon polarization index (summed over).")
    print(f"  v_n: The bosonic Matsubara frequency.")
    print(f"  {g_sq}: The square of the coupling constant 'g'.")
    print(f"  {q_j_sq}: The square of the j-th component of the momentum transfer vector q.")
    print(f"  {m}: The mass 'm' from the coupling term.")
    print(f"  {w_q_sq}: The square of the phonon frequency 'w_q' at momentum q.")

    print("-" * 60)

display_effective_interaction()