def reveal_hidden_word():
    """
    Solves the video game trivia and reveals the hidden word by combining
    the second letters of each answer.
    """
    answers = {
        1: "INVADE",
        2: "KONAMI",
        3: "EARTHWORM JIM",
        4: "GRIFFIN",
        5: "GAMING HEAD"
    }

    hidden_word = ""
    
    print("Revealing the hidden word by combining the second letters of each answer:")
    print("-" * 60)

    for i in range(1, 6):
        answer = answers[i]
        second_letter = answer[1]
        hidden_word += second_letter
        print(f"({i}) The answer is '{answer}'. The second letter is '{second_letter}'.")

    print("-" * 60)
    print(f"The hidden word is formed by the letters: { ' + '.join(list(hidden_word)) }")
    print(f"The final hidden word is: {hidden_word}")

reveal_hidden_word()