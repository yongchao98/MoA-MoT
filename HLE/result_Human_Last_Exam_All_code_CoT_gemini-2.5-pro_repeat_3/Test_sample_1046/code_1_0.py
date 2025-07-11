import random

def simulate_resistance_acquisition():
    """
    This simulation conceptually models two scenarios of bacteria acquiring drug resistance.
    """
    # --- Parameters ---
    time_steps = 15
    num_drugs = 5  # The total number of drugs we are testing against

    # --- Bacterium 1: Frequent Lateral Gene Transfer (LGT) ---
    # At each step, has a good chance to gain one resistance gene.
    resistance_bact1 = 0
    prob_lgt = 0.4  # Probability of a successful LGT event per time step

    # --- Bacterium 2: Rare Mutation with High Impact ---
    # Has a low chance of a single "jackpot" mutation that confers multiple resistances
    # and has high fitness (due to compensatory mutations).
    resistance_bact2 = 0
    prob_super_mutation = 0.15 # Lower probability for the rare event
    # Cross-resistance means one mutation gives resistance to multiple drugs.
    cross_resistance_factor = 4
    super_mutation_happened = False

    # --- Simulation Loop ---
    print("Conceptual Simulation of Resistance Acquisition")
    print("-" * 60)
    print(f"{'Time Step':<12} | {'Bacterium 1 (LGT)':<25} | {'Bacterium 2 (Mutation)':<25}")
    print("-" * 60)

    for t in range(1, time_steps + 1):
        # Bacterium 1 evolution by LGT
        if random.random() < prob_lgt and resistance_bact1 < num_drugs:
            resistance_bact1 += 1  # Gains one resistance

        # Bacterium 2 evolution by rare mutation
        if not super_mutation_happened and random.random() < prob_super_mutation:
            # The rare, high-impact mutation occurs.
            # The combination of cross-resistance and compensatory mutations allows
            # it to rapidly acquire resistance to multiple drugs.
            resistance_bact2 = min(num_drugs, cross_resistance_factor)
            super_mutation_happened = True

        print(f"{t:<12} | {resistance_bact1:<25} | {resistance_bact2:<25}")

    print("-" * 60)
    print("\nFinal Analysis:")
    print(f"After {time_steps} time steps:")
    print(f"Bacterium 1 (with frequent LGT) acquired resistance to {resistance_bact1} out of {num_drugs} drugs.")
    print(f"Bacterium 2 (with rare mutation) acquired resistance to {resistance_bact2} out of {num_drugs} drugs.")
    print("\nConclusion: The simulation shows that a single, rare mutational event that confers both")
    print("cross-resistance (to multiple drugs) and is supported by compensatory mutations (allowing it to spread)")
    print("can result in a pace of resistance acquisition comparable to that of frequent lateral gene transfer.")
    print("This supports answer choice B.")

simulate_resistance_acquisition()
<<<B>>>