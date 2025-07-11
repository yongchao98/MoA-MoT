def reveal_hidden_word():
    """
    Solves the video game trivia puzzle to reveal a hidden word.
    """
    # The answers to the five trivia questions
    answers = [
        "OCCUPY",       # (1) Second letter: C
        "KONAMI",       # (2) Second letter: O
        "EARTHWORM",    # (3) Second letter: A
        "MILES",        # (4) Second letter: I
        "TRAGIC WAND"   # (5) Second letter: R
    ]

    # Extract the second letter from each answer
    second_letters = [answer[1] for answer in answers]

    # The hidden word is the combination of these letters
    hidden_word = "".join(second_letters)

    # Display the process as an "equation" of letters
    equation = " + ".join(second_letters)
    print(f"Combining the second letters: {equation}")

    # Print the final revealed word
    print(f"The hidden word is: {hidden_word}")

reveal_hidden_word()