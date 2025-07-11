def reveal_hidden_word():
    """
    Solves a multi-part trivia puzzle to reveal a hidden word.
    The word is formed by the second letter of each answer, in order.
    """
    
    # Answers to the five trivia questions
    answer1 = "COLONIZES"
    answer2 = "SPORTS INTERACTIVE"
    answer3 = "EARTHWORM JIM"
    answer4 = "LEO"
    answer5 = "TRAGIC WAND"
    
    # Extract the second letter from each answer
    letter1 = answer1[1]
    letter2 = answer2[1]
    letter3 = answer3[1]
    letter4 = answer4[1]
    letter5 = answer5[1]
    
    # Combine the letters to form the hidden word
    hidden_word = letter1 + letter2 + letter3 + letter4 + letter5
    
    # Print the explanation and the final result
    print(f"The answer to trivia (1) is '{answer1}'. The second letter is '{letter1}'.")
    print(f"The answer to trivia (2) is '{answer2}'. The second letter is '{letter2}'.")
    print(f"The answer to trivia (3) is '{answer3}'. The second letter is '{letter3}'.")
    print(f"The answer to trivia (4) is '{answer4}'. The second letter is '{letter4}'.")
    print(f"The answer to trivia (5) is '{answer5}'. The second letter is '{letter5}'.")
    
    print("\nCombining the second letters in order gives the equation:")
    print(f"'{letter1}' + '{letter2}' + '{letter3}' + '{letter4}' + '{letter5}' = {hidden_word}")

reveal_hidden_word()
<<<OPAER>>>