def find_next_diagnostic_step():
    """
    Analyzes a clinical case and determines the best next diagnostic step.
    The clinical presentation strongly suggests allergic contact dermatitis from textiles.
    This is based on the history (new workout clothes) and the location of the rash
    (friction areas on axillary folds, sparing the vault where deodorant is applied).
    """

    question = "Which of the following is the best next step in diagnosis?"

    answer_choices = {
        'A': 'Skin biopsy',
        'B': 'KOH preparation',
        'C': 'Topical steroid',
        'D': 'Patch test',
        'E': 'None of the above'
    }

    correct_answer_key = 'D'
    correct_answer_value = answer_choices[correct_answer_key]

    print("Clinical Question:")
    print(question)
    print("\nAnswer Choices:")
    for key, value in answer_choices.items():
        print(f"{key}. {value}")

    print("\n---")
    print("Conclusion:")
    print(f"The best next step to confirm the diagnosis of allergic contact dermatitis and identify the specific allergen is a '{correct_answer_value}'.")
    print(f"Therefore, the correct option is {correct_answer_key}.")

find_next_diagnostic_step()