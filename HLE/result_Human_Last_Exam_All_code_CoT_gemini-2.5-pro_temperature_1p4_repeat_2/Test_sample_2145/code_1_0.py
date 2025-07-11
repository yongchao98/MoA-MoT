# Plan: Solve each trivia question and identify the first letter of the answer.
# Then, combine the letters to form the hidden word.

# (1) The movie with the answer sequence A, A, D, A, C, A, D, D, A is "Slumdog Millionaire".
answer1 = "Slumdog Millionaire"
letter1 = answer1[0]

# (2) In Soviet films, "Enemies" were to be shown on the left. In diagrams, they are often at the bottom (e.g., hell).
answer2 = "Enemies"
letter2 = answer2[0]

# (3) The TV series featuring a character who wants to eat cats (hence "Cheshire Salad", etc.) is "ALF".
answer3 = "ALF"
letter3 = answer3[0]

# (4) This is a pun-based riddle. "Karate" (the martial art) was banned in the USSR.
# "Showing" (a demonstration) of Karate was banned on a Karate mat. The "leaders" who
# came to power via "coups" (which also means "flips/turns" in Russian) are the karate masters.
answer4 = "Karate"
letter4 = answer4[0]

# The hidden word is formed by the first letter of each answer.
hidden_word = letter1 + letter2 + letter3 + letter4

# Print the components and the final result as an "equation".
print("The first letters of the answers are:")
print(f"(1) {answer1} -> {letter1}")
print(f"(2) {answer2} -> {letter2}")
print(f"(3) {answer3} -> {letter3}")
print(f"(4) {answer4} -> {letter4}")
print("\nCombining the letters reveals the hidden word:")
print(f"{letter1} + {letter2} + {letter3} + {letter4} = {hidden_word}")
