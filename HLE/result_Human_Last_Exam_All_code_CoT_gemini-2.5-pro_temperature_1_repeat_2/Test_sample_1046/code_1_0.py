import textwrap

def solve_bacterial_resistance_problem():
    """
    Analyzes the paradox of rapid drug resistance acquisition in a bacterium with a stable genome.
    """

    print("Analyzing the biological scenario:")
    print("-" * 35)

    scenario = {
        "Bacterium 1 (Fast Track)": "Has very common lateral transfer. Acquires resistance quickly.",
        "Bacterium 2 (The Paradox)": "Has a stable genome (no lateral transfer), but acquires resistance at an equal pace."
    }

    for bacterium, description in scenario.items():
        print(f"{bacterium}: {description}")

    print("\nCore Question: How can the mutation-driven pathway of Bacterium 2 be as fast as the lateral-transfer pathway of Bacterium 1?")
    
    explanation = """
    For a population relying on mutation to become resistant quickly, the initial mutation event is rare. To compensate, the resulting mutant strain must be exceptionally successful and spread through the population extremely rapidly once it emerges. We need to find the option that describes the creation of such a 'super strain'.
    """
    print(textwrap.fill(explanation, 80))

    print("\nEvaluating the necessary components for a 'super strain':")
    components = {
        "1. Resistance Mutation": "The initial, rare event that confers resistance.",
        "2. Compensatory Mutation": "A second mutation that counteracts any fitness cost from the first mutation, making the bacterium healthier and more competitive.",
        "3. Greatly Increased Fitness": "The combination of mutations makes the bacterium even more fit than the original, non-resistant strain.",
        "4. Cross-Resistance": "The mutations provide resistance to multiple drugs, a massive advantage."
    }
    
    for component, desc in components.items():
        print(f"- {component}: {desc}")

    conclusion = """
    Option B is the only choice that combines all these powerful evolutionary drivers. It describes an initial resistance mutation followed by compensatory mutations that not only restore but greatly increase fitness, and also lead to cross-resistance. This creates a highly advantageous strain that can sweep through the population, explaining the rapid pace of resistance acquisition.
    """
    print("\nConclusion:")
    print(textwrap.fill(conclusion, 80))

    final_answer_choice = "B"
    print(f"\nTherefore, the most comprehensive answer is B.")

solve_bacterial_resistance_problem()
<<<B>>>