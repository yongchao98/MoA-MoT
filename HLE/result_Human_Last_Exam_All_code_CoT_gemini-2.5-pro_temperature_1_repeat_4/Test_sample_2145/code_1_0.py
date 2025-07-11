def solve_trivia():
    """
    Solves a multi-part trivia question to find a hidden word.
    """
    # Answer to question (1)
    # The sequence of answers A, A, D, A C, A, D, D, A matches the one in the film "Slumdog Millionaire".
    answer1 = "Slumdog Millionaire"
    letter1 = answer1[0]

    # Answer to question (2)
    # In art history and film, Hell is depicted on the left and at the bottom of hierarchies.
    answer2 = "Hell"
    letter2 = answer2[0]

    # Answer to question (3)
    # The cat-themed dishes point to the cat-eating protagonist of the series "ALF".
    answer3 = "ALF"
    letter3 = answer3[0]

    # Answer to question (4)
    # The movie "The Dictator" was banned in several countries and its main character is a composite of various real-life dictators who came to power via coups.
    answer4 = "The Dictator"
    letter4 = answer4[0]

    hidden_word = letter1 + letter2 + letter3 + letter4

    print(f"(1) The answer is '{answer1}'. First letter: {letter1}")
    print(f"(2) The answer is '{answer2}'. First letter: {letter2}")
    print(f"(3) The answer is '{answer3}'. First letter: {letter3}")
    print(f"(4) The answer is '{answer4}'. First letter: {letter4}")
    print("\nCombining the first letters reveals the hidden word.")
    print(f"{letter1} + {letter2} + {letter3} + {letter4} = {hidden_word}")

solve_trivia()