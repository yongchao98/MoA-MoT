def solve_puzzle():
    """
    This function solves the trivia puzzle and prints the step-by-step
    derivation of the hidden word.
    """
    # Answers to the trivia questions
    answer1 = "Slumdog Millionaire"
    answer2 = "Enemies"
    answer3 = "Alf"
    answer4 = "Tintin"

    # Extract the first letter of each answer
    first_letter1 = answer1[0]
    first_letter2 = answer2[0]
    first_letter3 = answer3[0]
    first_letter4 = answer4[0]

    # Combine the letters to form the hidden word
    hidden_word = first_letter1 + first_letter2 + first_letter3 + first_letter4

    # Print the explanation and the final equation
    print(f"The answer to question (1) is '{answer1}'. The first letter is '{first_letter1}'.")
    print(f"The answer to question (2) is '{answer2}'. The first letter is '{first_letter2}'.")
    print(f"The answer to question (3) is '{answer3}'. The first letter is '{first_letter3}'.")
    print(f"The answer to question (4) is '{answer4}'. The first letter is '{first_letter4}'.")
    print("-" * 20)
    print("Combining the first letters forms the hidden word.")
    print(f"Final Equation: (1){first_letter1} + (2){first_letter2} + (3){first_letter3} + (4){first_letter4} = {hidden_word}")

solve_puzzle()