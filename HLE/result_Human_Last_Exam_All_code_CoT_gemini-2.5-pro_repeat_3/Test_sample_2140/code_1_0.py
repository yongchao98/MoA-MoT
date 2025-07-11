def reveal_hidden_word():
    """
    Solves the video game trivia to reveal a hidden word.
    The hidden word is formed by combining the second letter of each answer.
    """
    # The answers to the five trivia questions
    answers = [
        "invade",        # Q1: The Soviet Union invades the Moon.
        "SEGA",          # Q2: Apology came from the publisher, SEGA.
        "Earthworm Jim", # Q3: Game with the "Eat your vegetables" easter egg.
        "Leo",           # Q4: Latin for lion, companion to a "Panthera".
        "MAGIC WAND"     # Q5: Phrase describing Dennis Fong's mouse skill.
    ]

    hidden_word = ""
    print("Finding the hidden word by taking the second letter of each answer:")
    
    # Iterate through the answers and build the hidden word
    for i, answer in enumerate(answers):
        # Get the second letter (at index 1)
        second_letter = answer[1]
        hidden_word += second_letter
        print(f"Answer {i+1}: '{answer}' -> Second Letter: '{second_letter}'")

    print("\n------------------------------------")
    print(f"The hidden word is: {hidden_word}")
    print("------------------------------------")

reveal_hidden_word()
<<<neaea>>>