import textwrap

def explain_bacterial_resistance():
    """
    Analyzes the scenario of bacterial resistance acquisition and provides a reasoned explanation.
    """
    print("Analyzing the biological puzzle...")
    print("-" * 30)
    print("Scenario: Bacterium 1 (fast lateral transfer) and Bacterium 2 (stable genome, slow mutation) acquire resistance at an equal pace.\n")
    print("Objective: Find the most plausible biological mechanism for Bacterium 2's rapid adaptation.\n")

    # Step 1: Model the fitness of Bacterium 1
    print("Step 1: Modeling Bacterium 1 (Lateral Transfer)")
    print("This bacterium can quickly acquire a pre-made resistance gene. This is a very efficient process.")
    b1_base_fitness = 100
    b1_resistance_gain = 50
    b1_final_fitness = b1_base_fitness + b1_resistance_gain
    print(f"Hypothetical Fitness Equation for Bacterium 1: Base Fitness ({b1_base_fitness}) + Resistance Gene Gain ({b1_resistance_gain}) = {b1_final_fitness}\n")

    # Step 2: Model the fitness of Bacterium 2
    print("Step 2: Modeling Bacterium 2 (Mutation-driven)")
    print("This bacterium must rely on a rare internal mutation. This is a slow start.")
    print("However, to match the pace, subsequent events must provide a massive advantage.")
    b2_base_fitness = 100
    resistance_mutation_cost = -20  # Resistance often comes with a fitness cost
    
    # Let's consider the elements from Answer B
    compensatory_mutation_gain = 30  # A mutation that erases the cost and adds a bonus
    cross_resistance_bonus = 40      # The mutation also confers resistance to other drugs, a huge advantage
    
    print("The most complete explanation involves a cascade of events following the initial rare mutation:")
    print("  - An initial resistance mutation occurs (which may have a fitness cost).")
    print("  - A 'compensatory mutation' follows, not just removing the cost but GREATLY increasing overall fitness.")
    print("  - The same changes also lead to 'cross-resistance' (resistance to multiple drugs), which is a massive evolutionary shortcut.")
    
    b2_final_fitness = b2_base_fitness + resistance_mutation_cost + compensatory_mutation_gain + cross_resistance_bonus
    
    print("\nHypothetical Fitness Equation for Bacterium 2 based on option B:")
    # The final equation with each number as requested
    print(f"Final Fitness = Base Fitness ({b2_base_fitness}) + Initial Mutation Cost ({resistance_mutation_cost}) + Compensatory Gain ({compensatory_mutation_gain}) + Cross-Resistance Bonus ({cross_resistance_bonus})")
    print(f"Final Result = {b2_final_fitness}")

    # Step 3: Compare and conclude
    print("-" * 30)
    print("Conclusion:")
    print(f"Bacterium 1's final fitness is {b1_final_fitness}.")
    print(f"Bacterium 2's potential final fitness is {b2_final_fitness}.")
    
    conclusion_text = (
        "As the model shows, a combination of a powerful compensatory mutation "
        "and the bonus of cross-resistance can make a mutation-driven evolutionary "
        "pathway surprisingly explosive. Once the rare mutation occurs, the resulting "
        "strain can become so successful that it spreads rapidly, allowing it to "
        "match the pace of resistance acquisition seen in bacteria with high rates "
        "of lateral gene transfer. Therefore, option B provides the most comprehensive explanation."
    )
    
    print("\n" + textwrap.fill(conclusion_text, width=80))

explain_bacterial_resistance()
print("\n<<<B>>>")