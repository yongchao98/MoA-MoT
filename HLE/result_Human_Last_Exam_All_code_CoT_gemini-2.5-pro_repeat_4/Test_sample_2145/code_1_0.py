def solve_trivia_puzzle():
    """
    Solves a series of trivia questions and reveals a hidden word
    formed by the first letter of each answer.
    """
    answers = {
        1: "Slumdog Millionaire",
        2: "Sinners",
        3: "ALF",
        4: "Tintin"
    }

    hidden_word = ""
    print("Finding the hidden word by solving each trivia question:")
    print("-" * 55)

    for i in range(1, 5):
        answer = answers[i]
        first_letter = answer[0]
        hidden_word += first_letter
        print(f"({i}) The answer is '{answer}'. The first letter is '{first_letter}'.")

    print("-" * 55)
    print(f"The hidden word is revealed by combining the first letters: {hidden_word}")
    # The prompt asks to output each number in the final equation.
    # We will format the final output to show the composition of the word.
    final_equation = " + ".join(list(hidden_word))
    print(f"The final equation is: {final_equation} = {hidden_word}")

solve_trivia_puzzle()
<<<SSAT>>>