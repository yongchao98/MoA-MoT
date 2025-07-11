import math

def calculate_neutrino_asymmetry():
    """
    Calculates the neutrino-antineutrino asymmetry from neutral kaon decay.

    This function demonstrates how CP violation in the decay of the long-lived
    neutral kaon (K_L^0) leads to an unequal number of neutrinos and
    antineutrinos, even when starting from a symmetric state.

    This happens because the decay rate of K_L^0 into positrons and electron-neutrinos
    is slightly different from its decay rate into electrons and electron-antineutrinos.
    """

    # --- Physics Parameters (from Particle Data Group) ---

    # A hypothetical number of long-lived neutral kaons (K_L^0) available to decay.
    # We choose a large number to see a clear statistical difference.
    num_kaons = 1_000_000

    # Total Branching Ratio (BR) for semileptonic K_L^0 decays (K_L^0 -> pi e nu).
    # This is the fraction of kaons that will decay via the relevant channels.
    br_semileptonic = 0.4055

    # The charge asymmetry parameter (delta_L). This measured value quantifies
    # the extent of CP violation in these decays.
    # delta_L = (Rate(e+) - Rate(e-)) / (Rate(e+) + Rate(e-))
    delta_L = 0.00332

    # --- Calculation ---

    # 1. Total number of kaons that will undergo semileptonic decay.
    num_semileptonic_decays = num_kaons * br_semileptonic

    # 2. The fraction of these decays producing neutrinos (associated with positrons, e+).
    # This fraction is slightly larger than 0.5 due to CP violation.
    fraction_neutrinos = (1 + delta_L) / 2

    # 3. The fraction of these decays producing antineutrinos (associated with electrons, e-).
    # This fraction is slightly smaller than 0.5.
    fraction_antineutrinos = (1 - delta_L) / 2

    # 4. Calculate the absolute number of neutrinos and antineutrinos produced.
    num_neutrinos = num_semileptonic_decays * fraction_neutrinos
    num_antineutrinos = num_semileptonic_decays * fraction_antineutrinos
    net_asymmetry = num_neutrinos - num_antineutrinos

    # --- Output ---
    print("--- Neutrino Asymmetry from Kaon Decay ---")
    print("\nAnswering the question: Yes, a neutrino-antineutrino asymmetry can be induced.")
    print("This occurs because of a fundamental property called CP Violation in the weak force,")
    print("which causes the decay rates into neutrinos and antineutrinos to be slightly different.")
    print("\nHere is a calculation based on experimental data:\n")
    print(f"Starting with {num_kaons:,} long-lived neutral kaons (K_L^0):")
    print(f"Total undergoing semileptonic decay: {num_kaons} * {br_semileptonic} = {num_semileptonic_decays:,.0f}\n")

    # Final equation format as requested
    print("Number of electron neutrinos (v_e) produced:")
    print(f"{num_semileptonic_decays:,.0f} * (1 + {delta_L}) / 2 = {num_neutrinos:,.3f}")

    print("\nNumber of electron antineutrinos (v_e_bar) produced:")
    print(f"{num_semileptonic_decays:,.0f} * (1 - {delta_L}) / 2 = {num_antineutrinos:,.3f}")

    print("\n--- Result ---")
    print("The final count is asymmetric:")
    # Using math.floor to present integer particles for the final inequality
    print(f"Equation: {math.floor(num_neutrinos)} (v_e) != {math.floor(num_antineutrinos)} (v_e_bar)")
    print(f"This created a net excess of {net_asymmetry:,.3f} neutrinos.")


# Run the calculation and print the results
calculate_neutrino_asymmetry()
