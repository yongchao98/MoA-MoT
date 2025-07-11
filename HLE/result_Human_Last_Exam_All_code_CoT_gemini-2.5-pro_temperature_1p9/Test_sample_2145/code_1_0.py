def solve_trivia_puzzle():
    """
    Solves a four-part trivia puzzle to find a hidden word.
    The first letter of the answer to each question is taken consecutively
    to spell out the hidden word.
    """

    # The answers to the four trivia questions
    answers = {
        1: "Slumdog Millionaire",
        2: "Commoners",
        3: "Alf",
        4: "Donald Duck"
    }

    print("Solving the puzzle step-by-step:\n")

    first_letters = []
    for num, answer in answers.items():
        first_letter = answer[0]
        first_letters.append(first_letter)
        print(f"({num}) The answer is '{answer}'. The first letter is '{first_letter}'.")

    hidden_word = "".join(first_letters)
    equation = " + ".join(first_letters)

    print("\nAssembling the hidden word from the first letters:")
    print(f"{equation} = {hidden_word}")

    # Final answer in the specified format is handled outside this function call.
    return hidden_word

# Run the solver and print the final result
final_word = solve_trivia_puzzle()
print(f"\nThe hidden word is: {final_word}")
print(f"<<<{final_word}>>>")
