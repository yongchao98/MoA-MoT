def explain_chemical_potential_limit():
    """
    Explains the fundamental limit on the chemical potential for bosons
    in the context of Bose-Einstein condensation.
    """
    
    print("1. The starting point is the Bose-Einstein distribution for the occupation number <n_i> of an energy state epsilon_i:")
    print("   <n_i> = 1 / (exp[(epsilon_i - mu) / (k_B * T)] - 1)")
    print("   where 'mu' is the chemical potential.\n")

    print("2. For the occupation number <n_i> to be positive and non-infinite, the denominator must be positive.")
    print("   exp[(epsilon_i - mu) / (k_B * T)] - 1 > 0")
    print("   This simplifies to the condition: mu < epsilon_i for all energy states i.\n")

    print("3. The most stringent constraint is set by the lowest energy state, the ground state (epsilon_0).")
    print("   Therefore, the fundamental limit is: mu <= epsilon_0\n")

    print("4. During Bose-Einstein condensation (at temperatures T <= T_c), a macroscopic number of particles occupy the ground state.")
    print("   To accommodate this, the chemical potential 'mu' gets 'pinned' at the ground state energy, approaching it from below.")
    print("   In the idealized thermodynamic limit, mu = epsilon_0 for T <= T_c.\n")
    
    print("5. By definition, the chemical potential of a non-interacting Bose gas at absolute zero (T=0) is precisely the ground state energy, epsilon_0.")
    print("   At T=0, all particles are in the ground state, and adding one more costs epsilon_0 energy.\n")

    print("Conclusion: The fundamental limit is that the chemical potential 'mu' in a condensate must be equal to the ground state energy 'epsilon_0', which is the same as the chemical potential of a non-interacting Bose gas at zero temperature.\n")

    print("Therefore, the correct answer is C.")

explain_chemical_potential_limit()