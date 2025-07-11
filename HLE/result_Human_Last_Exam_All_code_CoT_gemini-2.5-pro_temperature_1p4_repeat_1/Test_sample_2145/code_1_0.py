def solve_trivia_puzzle():
    """
    This function solves the multi-part trivia puzzle and reveals the hidden word.
    """

    # First letter from the answer to question 1: (S)lumdog Millionaire
    letter_1 = 'S'

    # First letter from the answer to question 2: (E)nemies / (E)vil
    letter_2 = 'E'

    # First letter from the answer to question 3: (A)LF
    letter_3 = 'A'

    # First letter from the answer to question 4: (C)hipollino
    letter_4 = 'C'

    # The hidden word is formed by concatenating the first letters of each answer.
    hidden_word = letter_1 + letter_2 + letter_3 + letter_4

    print(f"The first letter of 'Slumdog Millionaire' is: {letter_1}")
    print(f"The first letter of 'Enemies' is: {letter_2}")
    print(f"The first letter of 'ALF' is: {letter_3}")
    print(f"The first letter of 'Chipollino' is: {letter_4}")
    print(f"The hidden word is: {hidden_word}")

solve_trivia_puzzle()