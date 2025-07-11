import random

def simulate_resistance_acquisition(generations=50):
    """
    Models two scenarios of acquiring drug resistance to illustrate the concept.
    """

    print("--- Scenario 1: Lateral Gene Transfer (LGT) ---")
    print("This bacterium acquires pre-existing resistance genes from others.")
    lgt_resistances = 0
    # LGT events happen at certain intervals
    lgt_event_generations = [10, 25, 40]
    for gen in lgt_event_generations:
        previous_resistances = lgt_resistances
        lgt_resistances += 1
        print(f"Generation {gen}: An LGT event occurred.")
        print(f"  Equation: Final Resistances ({lgt_resistances}) = Previous Resistances ({previous_resistances}) + New Gene (1)")

    print("\n--- Scenario 2: Mutation, Compensation, and Cross-Resistance ---")
    print("This bacterium relies on mutation within its stable genome.")
    mut_resistances = 0
    mut_fitness = 1.0

    # A rare resistance mutation occurs
    gen_resistance_mut = 12
    previous_resistances_mut = mut_resistances
    mut_resistances += 1
    mut_fitness *= 0.7 # This mutation comes with a 30% fitness cost
    print(f"Generation {gen_resistance_mut}: A rare resistance mutation occurred.")
    print(f"  Equation: Final Resistances ({mut_resistances}) = Previous Resistances ({previous_resistances_mut}) + New Mutation (1)")
    print(f"  Result: Fitness dropped to {mut_fitness:.2f} due to the cost of resistance.")


    # A compensatory mutation follows, which also provides cross-resistance
    gen_compensatory_mut = 18
    previous_resistances_mut = mut_resistances
    # This single mutation adds 2 new resistances (cross-resistance)
    cross_resistance_gain = 2
    newly_acquired = cross_resistance_gain
    mut_resistances += newly_acquired
    # It also restores and increases fitness
    mut_fitness = 1.2
    print(f"Generation {gen_compensatory_mut}: A compensatory mutation occurred.")
    print(f"  Effect: Gained cross-resistance to {newly_acquired} other drugs and fitness was restored/enhanced.")
    print(f"  Equation: Final Resistances ({mut_resistances}) = Previous Resistances ({previous_resistances_mut}) + Cross-Resistances ({newly_acquired})")
    print(f"  Result: Fitness increased to {mut_fitness:.2f}, allowing the strain to spread rapidly.")

    print("\n--- Conclusion ---")
    print(f"After ~{gen_compensatory_mut} generations, the mutation-based scenario yielded {mut_resistances} resistances.")
    print(f"The LGT scenario took {lgt_event_generations[1]} generations to achieve only 2 resistances.")
    print("This demonstrates that a highly advantageous mutational pathway (including compensation and cross-resistance) can indeed match the pace of LGT.")


simulate_resistance_acquisition()
