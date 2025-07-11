# 1. Define the answers to each trivia question.
answer1 = "Rocky"
answer2 = "Oppressed"
answer3 = "Alice in Borderland"
answer4 = "Nineteen Eighty-Four"

# 2. Extract the first letter from each answer.
letter1 = answer1[0]
letter2 = answer2[0]
letter3 = answer3[0]
letter4 = answer4[0]

# 3. Combine the letters to reveal the hidden word.
hidden_word = letter1 + letter2 + letter3 + letter4

# 4. Print the components of the "final equation" and the result.
print("The hidden word can be found by taking the first letter of each answer:")
print(f"(1) The answer is '{answer1}', and the first letter is '{letter1}'.")
print(f"(2) The answer is '{answer2}', and the first letter is '{letter2}'.")
print(f"(3) The answer is '{answer3}', and the first letter is '{letter3}'.")
print(f"(4) The answer is '{answer4}', and the first letter is '{letter4}'.")
print("\nCombining these letters forms the final equation:")
print(f"{letter1} + {letter2} + {letter3} + {letter4} = {hidden_word}")