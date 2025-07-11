def reveal_hidden_word():
    """
    Solves a series of trivia questions and combines the second letter
    of each answer to reveal a hidden word.
    """
    # The answers to the five trivia questions, derived from the logic above.
    answer1 = "CLAIMS"
    answer2 = "KONAMI"
    answer3 = "OGRE"
    answer4 = "MILES"
    answer5 = "PC GAMING"

    answers = [answer1, answer2, answer3, answer4, answer5]
    
    # Extract the second letter (index 1) from each answer.
    second_letters = [ans[1] for ans in answers]
    
    hidden_word = "".join(second_letters)

    # Print the equation showing how the word was formed.
    equation = " + ".join(second_letters) + " = " + hidden_word
    print(f"The answers to the trivia questions are: {', '.join(answers)}")
    print("Combining the second letter of each answer:")
    print(equation)
    print(f"\nThe hidden word is: {hidden_word}")

reveal_hidden_word()