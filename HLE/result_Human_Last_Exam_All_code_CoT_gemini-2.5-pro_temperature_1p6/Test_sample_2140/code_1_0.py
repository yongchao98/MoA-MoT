# The answers to the five trivia questions
answer1 = "INVADES"
answer2 = "EA"
answer3 = "EARTHWORM JIM"
answer4 = "VIR"
answer5 = "GAMING HAND"

# A list of all the answers
answers = [answer1, answer2, answer3, answer4, answer5]

# A list to store the second letters
second_letters = []

# Loop through each answer to get the second letter and print the process
for i, answer in enumerate(answers):
    # The second letter is at index 1
    second_letter = answer[1]
    second_letters.append(second_letter)
    print(f"The answer to question ({i+1}) is '{answer}'. The second letter is '{second_letter}'.")

# Combine the letters to form the hidden word
hidden_word = "".join(second_letters)

# Print the final result
print("\nCombining the second letters in order reveals the hidden word:")
print(hidden_word)