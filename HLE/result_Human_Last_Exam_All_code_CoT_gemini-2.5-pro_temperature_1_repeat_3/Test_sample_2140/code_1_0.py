# The answers to the five trivia questions
answer1 = "ESTABLISH"
answer2 = "KONAMI"
answer3 = "EARTHWORM JIM"
answer4 = "GRIFFIN"
answer5 = "ESDF"

# A list to hold all the answers
answers = [answer1, answer2, answer3, answer4, answer5]

# A list to hold the second letter of each answer
second_letters = []

# Loop through each answer and get the second letter (at index 1)
for answer in answers:
    if len(answer) > 1:
        second_letters.append(answer[1])

# Join the letters to form the hidden word
hidden_word = "".join(second_letters)

# Print the final result
print(f"The second letters are: {', '.join(second_letters)}")
print(f"The hidden word is: {hidden_word}")
