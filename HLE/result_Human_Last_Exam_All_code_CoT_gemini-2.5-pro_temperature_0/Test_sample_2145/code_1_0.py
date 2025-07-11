def solve_trivia_puzzle():
    """
    Solves a multi-part trivia puzzle to find a hidden word.
    """
    # Question 1
    answer1 = "Slumdog Millionaire"
    first_letter1 = answer1[0]
    print(f"(1) The Oscar-winning film is '{answer1}'. The first letter is {first_letter1}.")

    # Question 2
    answer2 = "Sinners"
    first_letter2 = answer2[0]
    print(f"(2) THEY are '{answer2}'. The first letter is {first_letter2}.")

    # Question 3
    answer3 = "ALF"
    first_letter3 = answer3[0]
    print(f"(3) The TV series is '{answer3}'. The first letter is {first_letter3}.")

    # Question 4
    answer4 = "Chilly Willy"
    first_letter4 = answer4[0]
    print(f"(4) X is the cartoon '{answer4}' (banned in Chile). The first letter is {first_letter4}.")

    # Reveal the hidden word
    hidden_word = first_letter1 + first_letter2 + first_letter3 + first_letter4
    print(f"\nThe hidden word is formed by the first letters: {first_letter1} + {first_letter2} + {first_letter3} + {first_letter4}")
    print(f"The hidden word is: {hidden_word}")

solve_trivia_puzzle()