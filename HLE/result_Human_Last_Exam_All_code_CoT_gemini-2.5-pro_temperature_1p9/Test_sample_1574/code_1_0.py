def solve_bulgakov_riddle():
    """
    Solves the literary riddle about Mikhail Bulgakov's "A Dog's Heart"
    and demonstrates the calculation mentioned in the prompt.
    """
    # Details from the prompt
    year_confiscated = 1926
    years_later_published = 40

    # Fulfilling the "equation" requirement
    year_published = year_confiscated + years_later_published
    print(f"Based on the prompt, the novel confiscated in {year_confiscated} was published {years_later_published} years later.")
    print(f"Equation: {year_confiscated} + {years_later_published} = {year_published}")
    print("-" * 20)

    # Answering the main question
    answer_choices = {
        "A": "Vasnetsova",
        "B": "Varvana",
        "C": "Zina",
        "D": "Maria",
        "E": "Darya"
    }

    correct_option = "A"
    correct_character_name = answer_choices[correct_option]

    # Explanation
    print("In the story, Polygraf Polygrafovich Sharikov brings home a typist, Comrade Vasnetsova.")
    print("He attempts to force her to marry him by threatening her with being fired from her job.")
    print("This act of coercion is the attempted assault in question.")
    print(f"\nThe correct answer is {correct_option}: {correct_character_name}.")

solve_bulgakov_riddle()