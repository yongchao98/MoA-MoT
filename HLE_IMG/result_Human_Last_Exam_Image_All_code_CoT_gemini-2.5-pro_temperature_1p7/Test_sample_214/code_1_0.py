import sys

# In Python, we can represent the lines of the poem and the syllable counts of each word.
# This approach is more reliable for a specific poem than a generic syllable-counting algorithm.
syllable_counts = {
    "ghostly": 2,
    "velum": 2,
    "forms": 1,
    "like": 1,
    "a": 1,
    "dance": 1,
    "vortex": 2,
    "nacreous": 3,
    "wavers": 2
}

# The actual poem is structured in three lines, based on the artist's original work.
line1_words = ["ghostly", "velum"]
line2_words = ["forms", "like", "a", "dance", "a", "vortex"]
line3_words = ["nacreous", "wavers"]

poem = [line1_words, line2_words, line3_words]
poem_lines_text = [" ".join(line) for line in poem]

print("Analyzing the syllable structure of the poem...")
print("-" * 40)

# Calculate and print the syllable count for each line
syllable_structure = []
for i, line_words in enumerate(poem):
    line_syllables = sum(syllable_counts[word] for word in line_words)
    syllable_structure.append(line_syllables)
    
    # Construct the equation string as requested
    equation_parts = [str(syllable_counts[word]) for word in line_words]
    equation_str = " + ".join(equation_parts)

    print(f"Line {i+1}: '{poem_lines_text[i]}'")
    print(f"Syllable equation: {equation_str} = {line_syllables}")
    print("-" * 40)

final_structure = "-".join(map(str, syllable_structure))
print(f"The poem follows a {final_structure} syllable structure.")
print("This is a common variation of the 5-7-5 syllable structure used in English-language Haiku.")
print("\nThe poetic form is Haiku.")
