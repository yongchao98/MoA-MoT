# This script determines the hidden word by solving each trivia question.

# Question 1: The film where a character's path to success is shown by the
# sequence of game show answers A, A, D, A, C, A, D, D, A.
# Answer: Slumdog Millionaire
letter1 = "S"

# Question 2: The group that Stalin would want on the left in films (political left)
# and that is at the bottom of historical diagrams (social pyramid).
# Answer: Lower Classes
letter2 = "L"

# Question 3: The TV series where cat-based dishes would be favorites.
# Answer: ALF
letter3 = "A"

# Question 4: The name 'X' for a film/series that was banned by 'X' (a play on words).
# Answer: Vlast (The Russian word for 'Power', a film banned by the 'powers that be').
letter4 = "V"

# Revealing the hidden word by consecutively naming the letters.
print(f"The letter for question (1) is: {letter1}")
print(f"The letter for question (2) is: {letter2}")
print(f"The letter for question (3) is: {letter3}")
print(f"The letter for question (4) is: {letter4}")

# The final hidden word.
hidden_word = letter1 + letter2 + letter3 + letter4
print(f"\nThe hidden word is: {hidden_word}")