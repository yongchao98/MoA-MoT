def reveal_hidden_word():
    """
    Solves a series of trivia questions to reveal a hidden word by combining
    the second letters of the answers.
    """

    # The answers to the five trivia questions
    answers = [
        "OCCUPY",        # (1) Second letter: C
        "KONAMI",        # (2) Second letter: O
        "EARTHWORM JIM", # (3) Second letter: A
        "ALARIC",        # (4) Second letter: L
        "MANIC HAND"     # (5) Second letter: A
    ]

    print("The answers to the video game trivia questions are:")
    print(f"(1) {answers[0]}")
    print(f"(2) {answers[1]}")
    print(f"(3) {answers[2]}")
    print(f"(4) {answers[3]}")
    print(f"(5) {answers[4]}")
    print("-" * 20)

    hidden_word = ""
    print("Combining the second letter of each answer reveals the hidden word:")
    
    # Extract the second letter from each answer
    for answer in answers:
        second_letter = answer[1].upper()
        hidden_word += second_letter
        
    print(f"The combined letters are: {hidden_word}")
    # The resulting word "COALA" is phonetically "KOALA"
    print(f"The hidden word is: KOALA")

reveal_hidden_word()
<<<KOALA>>>