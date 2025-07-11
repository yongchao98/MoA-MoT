def solve_trivia_puzzle():
    """
    Solves the multi-part trivia puzzle and reveals the hidden word.
    """

    # Answer to question 1
    answer1 = "Slumdog Millionaire"
    first_letter1 = answer1[0]
    print(f"(1) The film is '{answer1}'. The first letter is '{first_letter1}'.")

    # Answer to question 2
    answer2 = "Sinners"
    first_letter2 = answer2[0]
    print(f"(2) 'THEY' are '{answer2}'. The first letter is '{first_letter2}'.")

    # Answer to question 3
    answer3 = "ALF"
    first_letter3 = answer3[0]
    print(f"(3) The series is '{answer3}'. The first letter is '{first_letter3}'.")

    # Answer to question 4
    answer4 = "Kukly"
    first_letter4 = answer4[0]
    print(f"(4) 'X' is the show '{answer4}'. The first letter is '{first_letter4}'.")

    # The hidden word
    hidden_word = first_letter1 + first_letter2 + first_letter3 + first_letter4
    print(f"\nThe hidden word is formed by the first letters: {hidden_word}")

solve_trivia_puzzle()