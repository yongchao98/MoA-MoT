import math

def calculate_asymmetry():
    """
    Calculates and demonstrates the neutrino-antineutrino asymmetry from K_L decay.

    This is due to CP violation in the neutral kaon system, which causes the decay rate
    of K_L into neutrinos to be different from the rate into antineutrinos.
    """
    # The experimentally measured semileptonic charge asymmetry for K_L decays.
    # This parameter measures the fractional difference between the two decay rates.
    # delta_L = (Rate(K_L -> e+ nu) - Rate(K_L -> e- anti-nu)) / (Sum of Rates)
    delta_L = 3.32e-3

    # Let's consider a large, illustrative number of total semileptonic K_L decays.
    total_semileptonic_decays = 1_000_000

    # From the definition of delta_L, we can derive the relative rates:
    # Rate(nu) = (1 + delta_L) / 2
    # Rate(anti-nu) = (1 - delta_L) / 2
    # We can use these to find the number of neutrinos and antineutrinos produced.

    num_neutrinos = total_semileptonic_decays * (1 + delta_L) / 2
    num_antineutrinos = total_semileptonic_decays * (1 - delta_L) / 2

    print("This script quantifies the neutrino-antineutrino asymmetry from neutral kaon decays.")
    print(f"The known CP-violating charge asymmetry (δ_L) is {delta_L}.\n")
    print(f"Assuming {total_semileptonic_decays:,} total semileptonic K_L decays:\n")

    # Print the calculation for the number of neutrinos
    print("Number of neutrinos (νe) produced:")
    print(f"{total_semileptonic_decays} * (1 + {delta_L}) / 2 = {math.ceil(num_neutrinos):,}")

    # Print the calculation for the number of antineutrinos
    print("\nNumber of antineutrinos (ν̅e) produced:")
    print(f"{total_semileptonic_decays} * (1 - {delta_L}) / 2 = {math.floor(num_antineutrinos):,}")

    # Print the resulting difference
    excess_neutrinos = math.ceil(num_neutrinos) - math.floor(num_antineutrinos)
    print("\n------------------------------------------------------------")
    print(f"Result: An excess of {excess_neutrinos:,} neutrinos is generated.")
    print("This shows that kaon decays can indeed induce an asymmetry.")

if __name__ == "__main__":
    calculate_asymmetry()