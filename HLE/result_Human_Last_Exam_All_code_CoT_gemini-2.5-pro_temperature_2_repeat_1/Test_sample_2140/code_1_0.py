# The solutions to the trivia questions
answer1 = "Infested"
answer2 = "UEFA"
answer3 = "Earthworm Jim"
answer4 = "Leo"
answer5 = "Panic Hand"

# A list of the answers
answers = [answer1, answer2, answer3, answer4, answer5]

# Extract the second letter from each answer
second_letters = [ans[1] for ans in answers]

# The hidden word is the combination of these letters
hidden_word = "".join(second_letters)

# Print the final result in the specified format
print(f"The second letters are: '{second_letters[0]}', '{second_letters[1]}', '{second_letters[2]}', '{second_letters[3]}', '{second_letters[4]}'.")
print(f"The hidden word is: {hidden_word}")
