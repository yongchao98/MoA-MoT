def solve_philosophy_question():
    """
    This function identifies and explains the philosopher from the Lviv-Warsaw School
    who distinguished between 'symptoms' and 'signals'.
    """
    philosopher = "Kazimierz Ajdukiewicz"
    school = "Lviv-Warsaw School"
    concept = "hallmark (or mark)"
    components = {
        "symptom": "a natural sign connected by a causal link (e.g., smoke is a symptom of fire)",
        "signal": "a conventional sign based on agreement or rule (e.g., a green light is a signal to go)"
    }

    print(f"The philosopher from the {school} who decomposes the notion of a '{concept}' is {philosopher}.")
    print("\nHe distinguished between:")
    for component, description in components.items():
        print(f"- {component.capitalize()}: {description}")

solve_philosophy_question()