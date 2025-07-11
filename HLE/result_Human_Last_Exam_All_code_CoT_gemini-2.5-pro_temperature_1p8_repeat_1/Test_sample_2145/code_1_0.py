def solve_trivia_puzzle():
    """
    Solves a multi-part trivia puzzle to reveal a hidden word.
    """
    # Answers to each trivia question
    answer1 = "Amadeus"
    answer2 = "Devils"
    answer3 = "Alf"
    answer4 = "Donald Duck"

    # Get the first letter of each answer
    first_letter1 = answer1[0]
    first_letter2 = answer2[0]
    first_letter3 = answer3[0]
    first_letter4 = answer4[0]

    # Combine the letters to form the hidden word
    hidden_word = first_letter1 + first_letter2 + first_letter3 + first_letter4

    # Print the step-by-step reveal
    print(f"(1) The answer is the movie '{answer1}'. The first letter is {first_letter1}.")
    print(f"(2) The answer is '{answer2}'. The first letter is {first_letter2}.")
    print(f"(3) The answer is the series '{answer3}'. The first letter is {first_letter3}.")
    print(f"(4) The answer is '{answer4}'. The first letter is {first_letter4}.")

    print("\nConsecutively naming the first letter of each answer reveals the hidden word:")
    print(f"{first_letter1} + {first_letter2} + {first_letter3} + {first_letter4} = {hidden_word}")


if __name__ == "__main__":
    solve_trivia_puzzle()
