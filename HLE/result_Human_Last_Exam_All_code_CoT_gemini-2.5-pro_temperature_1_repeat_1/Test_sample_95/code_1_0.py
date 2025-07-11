def solve_riddle():
    """
    This script solves the riddle by constructing the answer from its components.
    The riddle's clues point to a famous European mountain range.
    - Visible from Milan, but not from smoky northern cities.
    - The answer is the ALPS.
    We will form this word using its ASCII character codes as the "equation".
    """

    # The ASCII codes for the letters A, L, P, S.
    num1 = 65
    num2 = 76
    num3 = 80
    num4 = 83

    # The "equation" is the sequence of numbers that represents the answer.
    print(f"The final equation is composed of the numbers: {num1}, {num2}, {num3}, {num4}")

    # Convert the numbers back to characters to reveal the name.
    char1 = chr(num1)
    char2 = chr(num2)
    char3 = chr(num3)
    char4 = chr(num4)

    # Combine the characters to form the final answer.
    answer = char1 + char2 + char3 + char4

    print(f"The name of THEM is: {answer}")

solve_riddle()