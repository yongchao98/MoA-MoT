# Plan:
# 1. Identify the source words from the riddles:
#    - "sleeveless garments that drape...": ROBES
#    - "experienced and trusted individuals...": COUNSELORS
# 2. Form new words using letters from the source words:
#    - From ROBES -> SOBER
#    - From COUNSELORS -> COUNSEL
# 3. Combine the new words to form the ship name and print the components of the final word equation.

# The first word is derived from "ROBES".
word1 = "SOBER"

# The second word is derived from "COUNSELORS".
word2 = "COUNSEL"

# The final ship name from the Culture series is the combination of the two words.
final_ship_name = f"{word1} {word2}"

print("The solution is derived in two parts:")
print(f"1. The word for the garment is 'ROBES'. Using its letters, we form the word: '{word1}'")
print(f"2. The word for the advisors is 'COUNSELORS'. Using its letters, we form the word: '{word2}'")
print("\nCombining these gives the final ship name in the following equation:")

# The final instruction asks to output each part of the "equation".
# Here, we print each component word and the final combined name.
print(f"'{word1}' + '{word2}' = '{final_ship_name}'")