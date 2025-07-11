import math

def calculate_asymmetry():
    """
    Calculates the neutrino-antineutrino asymmetry from neutral kaon decays.

    This function demonstrates how CP violation in the decay of the long-lived
    neutral kaon (K_L) leads to an imbalance in the production of neutrinos
    and antineutrinos, even if kaons and anti-kaons are initially created
    in equal numbers.
    """
    # Define the total number of K_L semileptonic decays to simulate.
    # Semileptonic decays are those that produce leptons (like electrons and neutrinos).
    total_semileptonic_decays = 1_000_000

    # Define the experimentally measured charge asymmetry for K_L decays.
    # This value represents the fractional difference in the rates of the two
    # key decay modes. It is a direct measure of CP violation.
    # charge_asymmetry = [Rate(K_L -> π⁻ e⁺ ν) - Rate(K_L -> π⁺ e⁻ ν̅)] / [Sum of Rates]
    charge_asymmetry = 0.00334

    # The charge asymmetry allows us to find the number of neutrinos (N_nu) and
    # antineutrinos (N_antinu) produced.
    # From the definition: (N_nu - N_antinu) / total_decays = charge_asymmetry
    # And by definition: N_nu + N_antinu = total_decays
    #
    # We can solve this system of two linear equations:
    # Adding the two equations: 2 * N_nu = total_decays * (1 + charge_asymmetry)
    # Subtracting the first from the second: 2 * N_antinu = total_decays * (1 - charge_asymmetry)

    num_neutrinos = 0.5 * total_semileptonic_decays * (1 + charge_asymmetry)
    num_antineutrinos = 0.5 * total_semileptonic_decays * (1 - charge_asymmetry)

    # We use math.floor to get integer particle counts
    num_neutrinos = math.floor(num_neutrinos)
    num_antineutrinos = math.floor(num_antineutrinos)
    
    # Calculate the net number of particles
    net_neutrino_excess = num_neutrinos - num_antineutrinos

    print("Demonstration of Neutrino Asymmetry from Neutral Kaon (K_L) Decay")
    print("=" * 65)
    print(f"Even if kaons and antikaons are produced equally, CP violation in their decay can create an asymmetry.")
    print("\nRelevant K_L Decay Modes:")
    # Using unicode for nice output
    print("1. K_L → π⁻ + e⁺ + νₑ   (produces a neutrino)")
    print("2. K_L → π⁺ + e⁻ + ν̅ₑ   (produces an antineutrino)")

    print(f"\nThe experimentally measured charge asymmetry is: {charge_asymmetry}")

    print(f"\nSimulating for a total of {total_semileptonic_decays:,} decays:")
    print("-" * 65)
    print(f"Number of decays producing neutrinos (νₑ): {num_neutrinos:,}")
    print(f"Number of decays producing antineutrinos (ν̅ₑ): {num_antineutrinos:,}")
    print("-" * 65)
    
    print("\nFinal Asymmetry Equation:")
    print(f"{num_neutrinos:,} (neutrinos) - {num_antineutrinos:,} (antineutrinos) = {net_neutrino_excess:,} (net neutrino excess)")
    print("=" * 65)
    print("This shows that a net excess of neutrinos is generated.")

if __name__ == '__main__':
    calculate_asymmetry()
