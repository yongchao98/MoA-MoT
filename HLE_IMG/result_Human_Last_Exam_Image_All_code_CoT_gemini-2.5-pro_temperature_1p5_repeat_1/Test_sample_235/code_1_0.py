import nltk
from nltk.corpus import cmudict

# This script will count the syllables in the poem to identify its metric pattern.
# We hypothesize that the poem includes the title and page number, aiming for 
# the 17 syllables of an "American Sentence".

# Ensure the required NLTK data is available
try:
    arpabet = cmudict.dict()
except LookupError:
    print("NLTK's CMU Pronouncing Dictionary not found. Downloading...")
    nltk.download('cmudict', quiet=True)
    arpabet = cmudict.dict()

def count_syllables(word):
    """
    Counts syllables in a single word. It handles hyphenated words and
    uses the CMU Pronouncing Dictionary for accuracy.
    """
    # Prepare the word for lookup
    word = word.lower().strip("',.")
    
    if not word:
        return 0
    
    # Handle hyphenated words by summing the syllables of their parts
    if '-' in word:
        return sum(count_syllables(part) for part in word.split('-'))

    # Handle possessives by checking the base word
    if "'" in word:
        word = word.split("'")[0]

    # Look up the word in the dictionary
    if word in arpabet:
        # The number of syllables is the number of phonemes with a stress marker (a digit)
        # We take the first pronunciation listed.
        pronunciations = arpabet[word]
        return len([phoneme for phoneme in pronunciations[0] if phoneme[-1].isdigit()])
    else:
        # Basic fallback for words not in the dictionary (e.g., proper nouns, rare words)
        # This is a simple heuristic and may not be perfectly accurate.
        vowels = "aeiouy"
        count = 0
        if word[0] in vowels:
            count += 1
        for index in range(1, len(word)):
            if word[index] in vowels and word[index-1] not in vowels:
                count += 1
        if word.endswith("e"):
            count -= 1
        if count == 0:
            count = 1
        return count

# List of words forming the complete poem, including title and page number.
# The number "35" is converted to "thirty-five".
poem_words = [
    "The", "first", "step", 
    "thirty-five",
    "rules", "and", "lines,", "an", "intricate", "spider's", "web", "work"
]

syllable_counts = [count_syllables(word) for word in poem_words]

print("Syllable count for each word in the poem:")
for word, count in zip(poem_words, syllable_counts):
    print(f"- '{word.strip(',')}': {count} syllable(s)")

# Display the final equation and sum
total_syllables = sum(syllable_counts)
equation_str = " + ".join(map(str, syllable_counts))

print("\nFinal syllable calculation:")
print(f"{equation_str} = {total_syllables}")

# Conclude based on the total syllable count
if total_syllables == 17:
    print("\nThe poem is a single sentence with a total of 17 syllables.")
    print("This corresponds to the 'American Sentence' form.")
else:
    print(f"\nThe poem has a total of {total_syllables} syllables, which does not fit a known strict pattern.")
