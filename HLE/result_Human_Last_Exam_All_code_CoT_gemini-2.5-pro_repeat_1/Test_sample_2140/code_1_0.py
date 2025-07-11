def solve_puzzle():
    """
    Solves the video game trivia puzzle to reveal a hidden word.
    """
    # List of answers to the trivia questions
    answers = {
        1: "MINED",
        2: "KONAMI",
        3: "EARTHWORM JIM",
        4: "DEUS",
        5: "STATIC HAND"
    }

    # Extract the second letter from each answer
    second_letters = [ans[1] for ans in answers.values()]
    
    hidden_word = "".join(second_letters)

    # Print the results
    print("The answers to the trivia questions are:")
    for i, answer in answers.items():
        print(f"({i}) {answer}")
    
    print("\nCombining the second letter of each answer reveals the hidden word.")
    
    # Show the equation
    equation_parts = []
    for answer in answers.values():
        # Highlight the second letter
        highlighted_answer = answer[0] + answer[1].upper() + answer[2:]
        equation_parts.append(f"{highlighted_answer}[1] = {answer[1]}")
    
    print(" + ".join(equation_parts))

    print(f"\nThe hidden word is: {hidden_word}")

solve_puzzle()
<<<IOAET>>>