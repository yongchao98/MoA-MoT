def calculate_neutrino_asymmetry():
    """
    Calculates the number of neutrinos and antineutrinos produced from K_L decays,
    demonstrating the asymmetry caused by CP violation.
    """
    # The charge asymmetry (delta_L) in semileptonic K_L decays is a measured physical constant.
    # It's defined as: (Rate(nu) - Rate(anti-nu)) / (Rate(nu) + Rate(anti-nu))
    # A positive value means more neutrinos are produced than antineutrinos.
    # The experimental value is approximately +0.00332.
    charge_asymmetry_delta_L = 0.00332

    # Let's assume we have a large number of K_L particles that will decay
    # through these semileptonic channels.
    total_semileptonic_decays = 1_000_000

    # From the definition of asymmetry, we can derive the number of resulting neutrinos
    # and antineutrinos.
    # Let N_nu be the number of neutrinos and N_antinu be the number of antineutrinos.
    # 1) N_nu + N_antinu = total_semileptonic_decays
    # 2) (N_nu - N_antinu) / (N_nu + N_antinu) = charge_asymmetry_delta_L
    # From (2), N_nu - N_antinu = charge_asymmetry_delta_L * total_semileptonic_decays
    # Solving these two linear equations:
    # N_nu = total_semileptonic_decays * (1 + charge_asymmetry_delta_L) / 2
    # N_antinu = total_semileptonic_decays * (1 - charge_asymmetry_delta_L) / 2

    # Calculate the number of neutrinos produced
    num_neutrinos = total_semileptonic_decays * (1 + charge_asymmetry_delta_L) / 2

    # Calculate the number of antineutrinos produced
    num_antineutrinos = total_semileptonic_decays * (1 - charge_asymmetry_delta_L) / 2
    
    # The final net lepton number generated
    net_lepton_number = num_neutrinos - num_antineutrinos

    print("--- Asymmetry Calculation from K_L Decay ---")
    print(f"Starting with {total_semileptonic_decays:,} semileptonic K_L decays.")
    print(f"Using known charge asymmetry (delta_L) = {charge_asymmetry_delta_L}")
    print("\nCalculating the number of neutrinos (N_nu)...")
    print(f"Equation: N_nu = Total Decays * (1 + delta_L) / 2")
    # Outputting the numbers in the final equation as requested
    print(f"Calculation: {num_neutrinos:,.0f} = {total_semileptonic_decays:,} * (1 + {charge_asymmetry_delta_L}) / 2")
    
    print("\nCalculating the number of antineutrinos (N_antinu)...")
    print(f"Equation: N_antinu = Total Decays * (1 - delta_L) / 2")
    # Outputting the numbers in the final equation as requested
    print(f"Calculation: {num_antineutrinos:,.0f} = {total_semileptonic_decays:,} * (1 - {charge_asymmetry_delta_L}) / 2")

    print("\n--- Results ---")
    print(f"Number of neutrinos produced: {int(num_neutrinos)}")
    print(f"Number of antineutrinos produced: {int(num_antineutrinos)}")
    print(f"Net excess of neutrinos (N_nu - N_antinu): {int(net_lepton_number)}")
    print("\nConclusion: A net positive number of neutrinos is produced, creating an asymmetry.")

calculate_neutrino_asymmetry()