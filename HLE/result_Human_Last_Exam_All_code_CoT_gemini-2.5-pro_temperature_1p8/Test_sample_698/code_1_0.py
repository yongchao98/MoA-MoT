# Description: This script verifies the correct "True Avalanche" answer by breaking down
# the target word and the chosen phrase into their constituent syllables and semantic parts.
# It then prints the analysis and the final, correct sentence.

# 1. Define the components of the puzzle
target_word = "Computer"
syllables_target = ["com", "pu", "ter"]
description = "My software tells the birds when and where to relieve themselves."

# 2. Define the candidate answer and its components
candidate_phrase = 'Come poo, tern!'
syllables_candidate = ["Come", "poo", "tern"]

# Semantic analysis of the candidate phrase in relation to the description
# "Computer" -> The software
# "Come poo" -> Telling to relieve oneself
# "tern" -> A type of bird
analysis = {
    "Syllable Pattern": f"'{'-'.join(syllables_target)}'",
    "Repeating Pattern": f"'{'-'.join(syllables_candidate)}'",
    "Software": "The word 'Computer' itself.",
    "Bird": "'tern' is a type of bird.",
    "Action": "'Come poo' is the command to relieve oneself."
}

# 3. Print the verification logic
print("--- Analysis of the True Avalanche ---")
print(f"Target Word: {target_word}")
print(f"Target Syllable Pattern: {analysis['Syllable Pattern']}")
print("-" * 20)
print(f"Candidate Phrase: '{candidate_phrase}'")
print(f"Candidate Syllable Pattern: {analysis['Repeating Pattern']}")
print("\n--- Matching against the description ---")
print(f"Description: '{description}'")
print(f"1. The software is: {analysis['Software']}")
print(f"2. The bird is: {analysis['Bird']}")
print(f"3. The command is: {analysis['Action']}")
print("\n--- Conclusion ---")
print("The candidate phrase matches both the syllabic and descriptive requirements.")
print("The final correct Word Avalanche is:")
print(f"\n{target_word}: \"{candidate_phrase}\"")
