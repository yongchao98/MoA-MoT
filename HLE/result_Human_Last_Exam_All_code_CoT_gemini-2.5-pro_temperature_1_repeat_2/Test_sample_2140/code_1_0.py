def solve_puzzle():
    """
    Solves the video game trivia puzzle to reveal a hidden word.
    """
    # The answers to the five trivia questions
    answer1 = "MINED"
    answer2 = "KONAMI"
    answer3 = "EARTHWORM JIM"
    answer4 = "ASINUS"
    answer5 = "MOUSE GAME"

    answers = [answer1, answer2, answer3, answer4, answer5]

    # Extract the second letter from each answer
    second_letters = [ans[1] for ans in answers]

    # The hidden word is formed by combining these letters
    hidden_word = "".join(second_letters)

    # Print the explanation and the final equation
    print("The answers to the trivia questions are:")
    print(f"(1) {answer1}")
    print(f"(2) {answer2}")
    print(f"(3) {answer3}")
    print(f"(4) {answer4}")
    print(f"(5) {answer5}")
    print("\nCombining the second letter of each answer reveals the hidden word.")
    
    # Building the final equation string
    equation_str = " + ".join(second_letters) + f" = {hidden_word}"
    
    print("\nThe final equation is:")
    print(equation_str)
    
    print(f"\nThe hidden word is: {hidden_word}")

solve_puzzle()