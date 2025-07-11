# Plan:
# 1. Store the answers to the five trivia questions in a list.
# 2. Extract the second letter from each answer.
# 3. Combine the letters to form the hidden word.
# 4. Print the final word.

# The answers to the five questions in order.
answers = [
    "INVADE",       # Q1: What the Soviet Union does with the Moon.
    "SEGA",         # Q2: Who apologized for the anthem mix-up.
    "EARTHWORMJIM", # Q3: Game with the "Eat your vegetables..." easter egg. (Spacing removed for consistency)
    "LEGATUS",      # Q4: Latin word for the companion.
    "MAINHAND"      # Q5: Original phrase for "MAGIC HAND".
]

# Initialize an empty string to build the hidden word.
hidden_word = ""

# Loop through each answer to get the second letter.
# In Python, string indices start at 0, so the second letter is at index 1.
for answer in answers:
    second_letter = answer[1]
    hidden_word += second_letter

# Print the final revealed word.
print("The hidden word is: " + hidden_word)
