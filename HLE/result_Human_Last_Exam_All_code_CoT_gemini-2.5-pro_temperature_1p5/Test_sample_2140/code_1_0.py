# Plan:
# 1. Define the list of answers derived from the trivia questions.
# 2. Create an empty list to hold the explanation for each step.
# 3. Initialize an empty string that will become the hidden word.
# 4. Iterate through each answer. For each one, take the second character (at index 1).
# 5. Add a descriptive sentence to our details list and append the character to the hidden word string.
# 6. After the loop, print all the descriptive sentences.
# 7. Finally, print the resulting hidden word.

def reveal_hidden_word():
    """
    Solves the puzzle by finding the answers to the trivia,
    extracting the second letter of each answer, and combining them.
    """
    answers = [
        "Occupy",
        "EA Sports",
        "Earthworm Jim",
        "Tigris",
        "Mouse Hand"
    ]

    hidden_word_letters = []
    
    print("Building the hidden word by taking the second letter of each answer:")
    
    for i, answer in enumerate(answers):
        # The second letter is at index 1 in the string
        second_letter = answer[1].upper()
        hidden_word_letters.append(second_letter)
        print(f"({i+1}) The second letter of '{answer}' is '{second_letter}'")

    final_word = "".join(hidden_word_letters)
    print(f"\nCombining the letters in order gives the hidden word: {final_word}")

reveal_hidden_word()