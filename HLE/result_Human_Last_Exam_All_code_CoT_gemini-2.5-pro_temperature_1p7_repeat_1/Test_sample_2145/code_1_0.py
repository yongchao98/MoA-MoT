# The answers to the four trivia questions
answer1 = "Slumdog Millionaire"
answer2 = "Devils"
answer3 = "ALF"
answer4 = "Pinocchio"

# Extracting the first letter from each answer
letter1 = answer1[0]
letter2 = answer2[0]
letter3 = answer3[0]
letter4 = answer4[0]

# Combining the letters to form the hidden word
hidden_word = letter1 + letter2 + letter3 + letter4

# Printing the process like an equation to show how the word was formed
print(f"The first letter of the answer to question (1) is: {letter1}")
print(f"The first letter of the answer to question (2) is: {letter2}")
print(f"The first letter of the answer to question (3) is: {letter3}")
print(f"The first letter of the answer to question (4) is: {letter4}")
print("\nCombining these letters reveals the hidden word:")
print(f"{letter1} + {letter2} + {letter3} + {letter4} = {hidden_word}")
