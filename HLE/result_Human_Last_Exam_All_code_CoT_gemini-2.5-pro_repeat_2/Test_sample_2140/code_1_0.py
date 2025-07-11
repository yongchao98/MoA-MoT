def reveal_hidden_word():
    """
    Solves a multi-part trivia puzzle to reveal a hidden word.
    The hidden word is formed by taking the second letter of each answer.
    """

    # The answers to the five trivia questions
    answer1 = "WEAPONIZES"
    answer2 = "KONAMI"
    answer3 = "EARTHWORM JIM"
    answer4 = "LEO"
    answer5 = "GAMING HAND"

    answers = [answer1, answer2, answer3, answer4, answer5]
    hidden_word_letters = []

    print("Revealing the hidden word by combining the second letters of each answer:\n")

    # Loop through the answers, get the second letter, and print the process
    for i, answer in enumerate(answers, 1):
        second_letter = answer[1]
        hidden_word_letters.append(second_letter.upper())
        print(f"({i}) The answer is '{answer}'. The second letter is '{second_letter.upper()}'.")

    # Combine the letters to form the final hidden word
    hidden_word = "".join(hidden_word_letters)

    print("\n----------------------------------")
    print(f"The final hidden word is: {hidden_word}")
    print("----------------------------------")

if __name__ == '__main__':
    reveal_hidden_word()