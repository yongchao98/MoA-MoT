# The answers to the four trivia questions
answer1 = "Slumdog Millionaire"
answer2 = "Enemies"
answer3 = "Alice"
answer4 = "Tintin"

# Extracting the first letter of each answer
letter1 = answer1[0]
letter2 = answer2[0]
letter3 = answer3[0]
letter4 = answer4[0]

# Combine the letters to form the hidden word
hidden_word = letter1 + letter2 + letter3 + letter4

# Print the step-by-step solution
print(f"(1) The film is '{answer1}'. The first letter is '{letter1}'.")
print(f"(2) The answer is '{answer2}'. The first letter is '{letter2}'.")
print(f"(3) The series is '{answer3}'. The first letter is '{letter3}'.")
print(f"(4) The answer is '{answer4}'. The first letter is '{letter4}'.")
print("-" * 20)

# Print the final equation to reveal the hidden word
print("The final equation reveals the hidden word:")
print(f"({1}){letter1} + ({2}){letter2} + ({3}){letter3} + ({4}){letter4} = {hidden_word}")