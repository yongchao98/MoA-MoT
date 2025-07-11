def reveal_hidden_word():
    """
    This function solves a multi-part trivia puzzle to reveal a hidden word.
    It determines the answer to each question, finds the second letter of each answer,
    and combines them to form the final word.
    """

    # Answer to question (1)
    answer1 = "OCCUPIES"
    letter1 = answer1[1]
    print(f"The answer to question (1) is '{answer1}'. The second letter is '{letter1}'.")

    # Answer to question (2)
    answer2 = "SPORTS INTERACTIVE"
    letter2 = answer2[1]
    print(f"The answer to question (2) is '{answer2}'. The second letter is '{letter2}'.")

    # Answer to question (3)
    answer3 = "EARTHWORM JIM"
    letter3 = answer3[1]
    print(f"The answer to question (3) is '{answer3}'. The second letter is '{letter3}'.")

    # Answer to question (4)
    answer4 = "LEO"
    letter4 = answer4[1]
    print(f"The answer to question (4) is '{answer4}'. The second letter is '{letter4}'.")

    # Answer to question (5)
    answer5 = "QUAKE HAND"
    letter5 = answer5[1]
    print(f"The answer to question (5) is '{answer5}'. The second letter is '{letter5}'.")

    # Combine the letters to form the hidden word
    hidden_word = letter1 + letter2 + letter3 + letter4 + letter5
    print("\nCombining the second letters in order...")
    print(f"The final hidden word is: {letter1} + {letter2} + {letter3} + {letter4} + {letter5} = {hidden_word}")

reveal_hidden_word()
<<<CPAEU>>>