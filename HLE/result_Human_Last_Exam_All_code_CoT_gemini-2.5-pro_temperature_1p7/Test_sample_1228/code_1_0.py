def calculate_neutrino_asymmetry():
    """
    Calculates the neutrino-antineutrino asymmetry from the decay of
    long-lived neutral kaons (K_L), based on experimental data.
    """
    
    # --- Input Parameters ---
    
    # Let's start with a large, hypothetical number of K_L mesons that decay.
    # The actual number isn't important; this is for demonstrating the principle.
    initial_KL_mesons = 100_000_000
    
    # --- Experimental Data (from Particle Data Group) ---
    
    # Total branching ratio of K_L into semileptonic final states (e.g., π⁺e⁻ν, π⁻e⁺ν, etc.).
    # This combines decays to electrons/positrons (40.55%) and muons/antimuons (27.04%).
    br_semileptonic = 0.4055 + 0.2704
    
    # The charge asymmetry parameter (δ_L). This measures how much more often
    # K_L decays to positive leptons (l⁺ and ν) than negative leptons (l⁻ and ν̄).
    # δ_L = (Rate(ν) - Rate(ν̄)) / (Rate(ν) + Rate(ν̄))
    delta_L = 3.32e-3  # or 0.00332
    
    # --- Calculation ---
    
    # 1. Total number of K_L mesons that decay semileptonically
    num_semileptonic_decays = initial_KL_mesons * br_semileptonic
    
    # 2. We can solve for the number of neutrino and antineutrino producing decays.
    # Let N_ν be the number of decays producing a neutrino (K_L → l⁺ + ν + X).
    # Let N_ν̄ be the number of decays producing an antineutrino (K_L → l⁻ + ν̄ + X).
    # We have two equations:
    #   I) N_ν + N_ν̄ = num_semileptonic_decays
    #  II) (N_ν - N_ν̄) / (N_ν + N_ν̄) = delta_L  => N_ν - N_ν̄ = delta_L * num_semileptonic_decays
    # Solving this system gives:
    num_neutrinos_produced = num_semileptonic_decays * (1 + delta_L) / 2
    num_antineutrinos_produced = num_semileptonic_decays * (1 - delta_L) / 2
    
    # The net asymmetry is the difference between them.
    net_neutrino_number = num_neutrinos_produced - num_antineutrinos_produced
    
    # --- Output Results ---
    
    print(f"Starting with {int(initial_KL_mesons):,} K_L mesons.")
    print("-" * 50)
    print(f"Total semileptonic decays: {int(num_semileptonic_decays):,}")
    print(f"Charge asymmetry parameter δ_L used: {delta_L}")
    print("-" * 50)
    print("Final numbers from the decay process:")
    print(f"Number of neutrinos (ν) produced: {int(round(num_neutrinos_produced)):,}")
    print(f"Number of antineutrinos (ν̄) produced: {int(round(num_antineutrinos_produced)):,}")
    print("\n--- Final Asymmetry Equation ---")
    print(f"{int(round(num_neutrinos_produced))} neutrinos - {int(round(num_antineutrinos_produced))} antineutrinos = {int(round(net_neutrino_number)):,} net neutrinos")
    
if __name__ == "__main__":
    calculate_neutrino_asymmetry()
