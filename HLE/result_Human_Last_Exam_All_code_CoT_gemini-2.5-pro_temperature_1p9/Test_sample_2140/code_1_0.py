# The answers to the five trivia questions
answer1 = "BLOW UP"
answer2 = "SPORTS INTERACTIVE"
answer3 = "EARTHWORM JIM"
answer4 = "LEO"
answer5 = "TRAGIC SCAN"

# Store answers in a list
answers = [answer1, answer2, answer3, answer4, answer5]

# Extract the second letter from each answer
second_letters = [ans[1] for ans in answers]

# Build the hidden word
hidden_word = "".join(second_letters)

# Print the process and the final result
print("Revealing the hidden word by combining the second letters:")
print(f"1. {answers[0][:1]}[{answers[0][1]}]{answers[0][2:]}")
print(f"2. {answers[1][:1]}[{answers[1][1]}]{answers[1][2:]}")
print(f"3. {answers[2][:1]}[{answers[2][1]}]{answers[2][2:]}")
print(f"4. {answers[3][:1]}[{answers[3][1]}]{answers[3][2:]}")
print(f"5. {answers[4][:1]}[{answers[4][1]}]{answers[4][2:]}")
print("\nThe final hidden word is:")
print(hidden_word)
