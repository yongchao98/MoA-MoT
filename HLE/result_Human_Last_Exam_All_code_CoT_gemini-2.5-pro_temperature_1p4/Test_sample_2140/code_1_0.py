def reveal_hidden_word():
    """
    This script solves a multi-part trivia puzzle to reveal a hidden word.
    It identifies the answer to each trivia question and then combines the
    second letter of each answer to form the final word.
    """

    # The answers to the five trivia questions
    answer1 = "MINED"
    answer2 = "KONAMI"
    answer3 = "WARCRAFT"
    answer4 = "LEO"
    answer5 = "ROCKETJUMP"

    answers = [answer1, answer2, answer3, answer4, answer5]

    # Extract the second letter from each answer
    hidden_word_letters = [word[1] for word in answers]
    hidden_word = "".join(hidden_word_letters)

    # Prepare the equation for printing
    equation_parts = []
    for i, answer in enumerate(answers):
        letter = answer[1]
        # Formatting each part of the equation, e.g., "(1) MINED -> 'I'"
        equation_parts.append(f"({i+1}) {answer} -> '{letter}'")

    # Print the explanation and the final result
    print("The hidden word is revealed by taking the second letter of each answer:")
    print("\n".join(equation_parts))
    print("\nCombining the letters results in the final equation:")
    # Formatting the final equation, e.g., 'I' + 'O' + 'A' + 'E' + 'O' = IOAEO
    final_equation = " + ".join([f"'{l}'" for l in hidden_word_letters])
    print(f"{final_equation} = {hidden_word}")


reveal_hidden_word()
