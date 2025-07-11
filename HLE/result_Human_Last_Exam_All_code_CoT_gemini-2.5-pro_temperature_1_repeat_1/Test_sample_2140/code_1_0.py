def reveal_hidden_word():
    """
    This function reveals the hidden word by solving five video game trivia questions
    and combining the second letter of each answer.
    """

    # Answers to the trivia questions
    answer1 = "PAINT"
    answer2 = "KONAMI"
    answer3 = "DOOM"
    answer4 = "LEO"
    answer5 = "MAIN HAND"

    # Extract the second letter from each answer
    letter1 = answer1[1]
    letter2 = answer2[1]
    letter3 = answer3[1]
    letter4 = answer4[1]
    letter5 = answer5[1] # The second letter of the first word

    # Combine the letters to form the hidden word
    hidden_word = letter1 + letter2 + letter3 + letter4 + letter5

    # Print the explanation and the result
    print("Finding the hidden word by combining the second letters of each answer:")
    print(f"(1) The word is {answer1}. The second letter is {letter1}.")
    print(f"(2) The word is {answer2}. The second letter is {letter2}.")
    print(f"(3) The word is {answer3}. The second letter is {letter3}.")
    print(f"(4) The word is {answer4}. The second letter is {letter4}.")
    print(f"(5) The phrase is '{answer5}'. The second letter is {letter5}.")
    print("\nCombining these letters in order reveals the hidden word:")
    print(f"'{letter1}' + '{letter2}' + '{letter3}' + '{letter4}' + '{letter5}' = {hidden_word}")

reveal_hidden_word()
<<<AOOEA>>>