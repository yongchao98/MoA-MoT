def solve_hidden_word_puzzle():
    """
    Solves the video game trivia puzzle to reveal a hidden word.
    """

    # The answers to the five trivia questions
    answers = [
        "INFILTRATE",
        "EA SPORTS",
        "EARTHWORM JIM",
        "LEO",
        "MAJIC WAND"
    ]

    hidden_word = ""

    print("Revealing the hidden word by finding the second letter of each answer:")
    print("-" * 70)

    # Loop through each answer, extract the second letter, and build the hidden word
    for i, answer in enumerate(answers):
        # The index for the second letter is 1
        second_letter = answer[1]
        hidden_word += second_letter
        # The prompt requires printing each part of the "equation"
        print(f"({i+1}) The answer is '{answer}'. The second letter is '{second_letter}'.")

    print("-" * 70)
    print("The final hidden word is:")
    print(hidden_word)


solve_hidden_word_puzzle()