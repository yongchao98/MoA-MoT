# This script constructs and prints the fourteen-syllable phrase as requested.

# The final phrase describes abnormal cell growth in the Luys' body.
# It is designed to match the rhyme and 14-syllable scansion of
# a doubled "Bibbidi-Bobbidi-Boo" from Disney's Cinderella.
phrase = "Cells Are Multiplying Now Inside The Body Of Luys"

# Per the instructions, we output each number in the syllable "equation".
# Each number represents the syllable count of the corresponding word in the phrase.
syllable_equation = "1 + 1 + 4 + 1 + 2 + 1 + 2 + 1 + 1 = 14"

print("Syllable Breakdown:")
print(syllable_equation)
print("-" * 25)
print("Final Phrase:")

# Capitalize each word and print the result.
print(phrase.upper())