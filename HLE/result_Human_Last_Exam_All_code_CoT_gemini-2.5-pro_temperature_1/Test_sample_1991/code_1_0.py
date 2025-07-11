# This script will generate and display the requested phrase.

# The phrase is designed to be fourteen syllables, describing the specified condition,
# and matching the rhyme and scansion of a famous Disney song title.
final_phrase = "The Malignant Growth In Luysi's Body Is Atrocious"

# Print the final phrase.
print(final_phrase)

# The following lines fulfill the requirement to "output each number in the final equation".
# This is interpreted as showing the syllable count for each word, which sums to 14.
# The words and their syllable counts are:
# The (1), Malignant (3), Growth (1), In (1), Luysi's (2), Body (2), Is (1), Atrocious (3)
syllable_numbers = [1, 3, 1, 1, 2, 2, 1, 3]
equation_string = " + ".join(map(str, syllable_numbers))
total_syllables = sum(syllable_numbers)

# Print the final equation representing the syllable count.
print(f"Syllable Equation: {equation_string} = {total_syllables}")