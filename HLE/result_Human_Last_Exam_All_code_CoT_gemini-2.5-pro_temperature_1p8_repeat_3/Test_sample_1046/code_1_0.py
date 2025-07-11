def evaluate_resistance_scenarios():
    """
    This script models the pace of drug resistance acquisition in different scenarios
    to determine the most plausible explanation for the observed phenomenon.
    """

    print("Analyzing how a bacterium with a stable genome could acquire drug resistance as fast as one with lateral transfer.\n")

    # --- Baseline Parameters (Arbitrary units for illustration) ---
    # Pace of resistance acquisition via Horizontal Gene Transfer (HGT). This is our benchmark.
    hgt_pace = 90

    # Pace of resistance from a single, rare mutation spreading vertically.
    rare_mutation_pace = 15

    # Factors affecting the pace of resistance acquisition through vertical transfer.
    # Resistance mutations often have a cost, reducing the bacterium's fitness and slowing its spread.
    fitness_cost_penalty = -20
    # Compensatory mutations fix the fitness cost and can even provide a significant boost.
    strong_compensation_boost = 60
    # Basic compensation just restores fitness to normal.
    basic_compensation_boost = 20
    # Cross-resistance provides resistance to multiple drugs, a significant advantage.
    cross_resistance_boost = 35

    print(f"Benchmark Pace (Bacterium 1 with HGT): {hgt_pace}\n")
    print("Evaluating scenarios for Bacterium 2 (no HGT):\n")

    # Scenario A: Rare mutations only.
    pace_a = rare_mutation_pace
    print(f"Scenario A (Rare Mutations):")
    print(f"Calculation: Base Pace = {pace_a}")
    print(f"Resulting Pace: {pace_a}. This is much slower than HGT pace of {hgt_pace}.\n")

    # Scenario B: Rare mutations, strong compensation, and cross-resistance.
    pace_b = rare_mutation_pace + fitness_cost_penalty + strong_compensation_boost + cross_resistance_boost
    print(f"Scenario B (Rare Mutations + Strong Compensation + Cross-Resistance):")
    print(f"Calculation: {rare_mutation_pace} (Base) + {fitness_cost_penalty} (Fitness Cost) + {strong_compensation_boost} (Strong Compensation) + {cross_resistance_boost} (Cross-Resistance)")
    print(f"Resulting Pace: {pace_b}. This pace is comparable to the HGT pace of {hgt_pace}.\n")

    # Scenario C: Contamination. This is a procedural error, not a biological mechanism.
    print(f"Scenario C (Contamination):")
    print("This scenario invalidates the experiment's premise. Pace cannot be calculated.\n")

    # Scenario D: Rare mutations, cross-resistance, but NO compensation. The fitness cost remains.
    pace_d = rare_mutation_pace + fitness_cost_penalty + cross_resistance_boost
    print(f"Scenario D (Rare Mutations + Cross-Resistance, No Compensation):")
    print(f"Calculation: {rare_mutation_pace} (Base) + {fitness_cost_penalty} (Fitness Cost) + {cross_resistance_boost} (Cross-Resistance)")
    print(f"Resulting Pace: {pace_d}. The uncompensated fitness cost makes this pace too low.\n")

    # Scenario E: Rare mutations with basic compensation, but NO cross-resistance.
    pace_e = rare_mutation_pace + fitness_cost_penalty + basic_compensation_boost
    print(f"Scenario E (Rare Mutations + Basic Compensation):")
    print(f"Calculation: {rare_mutation_pace} (Base) + {fitness_cost_penalty} (Fitness Cost) + {basic_compensation_boost} (Basic Compensation)")
    print(f"Resulting Pace: {pace_e}. Without the cross-resistance bonus and strong compensation, this pace is too low.\n")

    print("--- Conclusion ---")
    print("The only scenario that produces a resistance acquisition pace comparable to HGT is Scenario B.")
    print("This is because the combination of a greatly increased fitness from compensatory mutations and the bonus of cross-resistance creates exceptionally strong positive selection.")

if __name__ == '__main__':
    evaluate_resistance_scenarios()