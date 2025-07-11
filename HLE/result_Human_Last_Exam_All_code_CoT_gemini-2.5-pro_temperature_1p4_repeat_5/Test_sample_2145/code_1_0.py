def solve_trivia_puzzle():
    """
    Solves a multi-part trivia puzzle to reveal a hidden word.
    The final output explains how the hidden word is formed from the first letter of each answer.
    """

    # (1) The first trivia answer
    answer1 = "Slumdog Millionaire"
    letter1 = answer1[0]
    print(f"(1) The sequence A, A, D, A, C, A, D, D, A represents the correct multiple-choice answers Jamal gives on his path to success in the Oscar-winning film '{answer1}'. The first letter is '{letter1}'.")

    # (2) The second trivia answer
    answer2 = "Underclass"
    letter2 = answer2[0]
    print(f"(2) In Soviet ideology, the working class (proletariat) was politically on the 'left'. In social hierarchy diagrams, this group forms the base at the bottom. The answer is '{answer2}'. The first letter is '{letter2}'.")

    # (3) The third trivia answer
    answer3 = "Red Dwarf"
    letter3 = answer3[0]
    print(f"(3) The cat-themed dishes point to the TV series '{answer3}', which features the character Cat. The 'missing letter' is a pun on the show's running gag about missing (lost) mail. The first letter is '{letter3}'.")

    # (4) The fourth trivia answer
    answer4 = "Firefox"
    letter4 = answer4[0]
    print(f"(4) 'X' is the 1982 movie '{answer4}'. This Cold War film, critical of the Soviet Union (a state born from a coup), was banned there. The question is a complex pun. The first letter is '{letter4}'.")

    print("\n---")
    print("Combining the first letter of each answer reveals the hidden word:")
    print(f"{letter1} + {letter2} + {letter3} + {letter4} = {letter1}{letter2}{letter3}{letter4}")
    print("---\n")

# Execute the function to print the solution
solve_trivia_puzzle()

# The final hidden word
print("<<<SURF>>>")