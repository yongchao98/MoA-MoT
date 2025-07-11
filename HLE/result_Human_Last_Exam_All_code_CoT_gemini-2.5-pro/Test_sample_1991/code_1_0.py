# The fourteen-syllable phrase describing the condition
phrase = "Apoptosis Failing In The Luys Is So Atrocious"
print(phrase)

# The "equation" representing the syllable count for each word in the phrase
# as per the user's special instruction.
words = phrase.split()
syllable_counts = [4, 2, 1, 1, 1, 1, 1, 3] # Syllables for each word
total_syllables = sum(syllable_counts)

# Build and print the equation string
equation_str = " + ".join(map(str, syllable_counts)) + f" = {total_syllables}"
print("Syllable Equation: " + equation_str)