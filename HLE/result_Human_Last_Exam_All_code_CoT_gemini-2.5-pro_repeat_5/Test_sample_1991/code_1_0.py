# Plan:
# 1. Define the words of the phrase and their corresponding syllable counts.
# 2. Create and print an "equation" to verify that the syllables sum to 14.
# 3. Construct the final phrase with each word capitalized.
# 4. Print the resulting phrase.

# A dictionary mapping each word of the phrase to its syllable count.
syllable_data = {
    "Luysian": 3,
    "Body": 2,
    "Cells": 1,
    "Proliferating": 5,
    "Without": 2,
    "End": 1
}

# Extract the words and their syllable counts.
words = list(syllable_data.keys())
syllable_counts = list(syllable_data.values())

# Create the equation string from the syllable counts.
equation_string = " + ".join(map(str, syllable_counts))
total_syllables = sum(syllable_counts)

# Print the syllable count verification as an equation.
print("Syllable Count Equation:")
print(f"{equation_string} = {total_syllables}")
print("-" * 25)

# Construct the final phrase by joining the capitalized words.
final_phrase = " ".join(words)

# Print the final phrase.
print("Final Phrase:")
print(final_phrase)