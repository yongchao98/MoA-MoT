import math

def main():
    """
    This script models the generation of a neutrino-antineutrino asymmetry
    from the CP-violating decay of long-lived neutral kaons (K_L).
    This is the crucial final step in the process described.
    """

    # --- Input Parameters ---
    # Let's assume the process creates a population of 1,000,000 K_L mesons
    # that are available to decay.
    N_KL = 1_000_000

    # --- Physical Constants ---
    # The experimentally measured charge asymmetry for semileptonic K_L decays.
    # This value represents the fractional difference in decay rates.
    # delta_L = (Rate(KL -> e+ nu) - Rate(KL -> e- nubar)) / (Rate(KL -> e+ nu) + Rate(KL -> e- nubar))
    delta_L = 3.32e-3

    # Branching ratio for K_L decaying into the electron+pion+neutrino channel.
    # BR_semileptonic_electron = BR(K_L -> pi+ e- nubar_e) + BR(K_L -> pi- e+ nu_e)
    BR_semileptonic_electron = 0.4055

    # --- Calculation ---
    # Total number of K_L mesons that will decay via this specific channel.
    num_semileptonic_decays = N_KL * BR_semileptonic_electron

    # The net difference (the asymmetry) is the total number of decays
    # multiplied by the asymmetry parameter.
    net_neutrino_number = num_semileptonic_decays * delta_L

    # From the definition of delta_L, we can solve for the individual numbers of
    # neutrinos and antineutrinos produced from the total decays.
    # Let N_v be neutrinos and N_vb be antineutrinos.
    # N_v + N_vb = num_semileptonic_decays
    # N_v - N_vb = net_neutrino_number
    # Solving this system of equations gives:
    num_neutrinos = (num_semileptonic_decays + net_neutrino_number) / 2
    num_antineutrinos = (num_semileptonic_decays - net_neutrino_number) / 2

    # --- Output ---
    print(f"Starting with a population of {N_KL:,} K_L mesons:")
    print(f"Total relevant decays (K_L -> pi e nu/nubar): {math.trunc(num_semileptonic_decays)}")
    print("\n--- Asymmetry Calculation ---")
    print("The CP-violating decay creates an unequal number of neutrinos and antineutrinos:")
    print(f"Number of neutrinos (ν_e) produced = {math.trunc(num_neutrinos)}")
    print(f"Number of antineutrinos (ν̅_e) produced = {math.trunc(num_antineutrinos)}")
    print("\n--- Final Equation for the Net Neutrino Asymmetry ---")
    # We display the final relationship using the calculated numbers.
    # Net Asymmetry = Number of Neutrinos - Number of Antineutrinos
    print(
        f"{math.trunc(net_neutrino_number)} = {math.trunc(num_neutrinos)} - {math.trunc(num_antineutrinos)}"
    )

if __name__ == "__main__":
    main()