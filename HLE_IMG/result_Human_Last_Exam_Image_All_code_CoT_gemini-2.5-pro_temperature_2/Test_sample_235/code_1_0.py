import nltk

def count_syllables(word_list):
    """Counts syllables for a list of words using the CMU Pronouncing Dictionary."""
    try:
        d = nltk.corpus.cmudict.dict()
    except LookupError:
        print("Downloading CMU Pronouncing Dictionary...")
        nltk.download('cmudict', quiet=True)
        d = nltk.corpus.cmudict.dict()

    total_syllables = 0
    equation_parts = []

    for word in word_list:
        clean_word = word.lower().strip(".,?!")
        # Handle specific cases and fallbacks
        if clean_word == "pi":
            count = 1
        elif clean_word == "spider's":
             # CMU dict doesn't have the possessive, so we use the base word.
            count = [len(list(y for y in x if y[-1].isdigit())) for x in d['spider']][0]
        elif clean_word not in d:
            print(f"Warning: '{word}' not found in dictionary. Syllable count may be estimated.")
            count = 0 # Should not happen with the given words
        else:
            # Count the vowels (marked by digits) in the pronunciation
            count = [len(list(y for y in x if y[-1].isdigit())) for x in d[clean_word]][0]

        total_syllables += count
        equation_parts.append(str(count))
        print(f"The word '{word}' has {count} syllable(s).")

    return equation_parts, total_syllables

# The complete text combined from the page's elements
# "The First Step Thirty-Five Pi rules and lines an intricate spider's web work"
poem_words = [
    "The", "First", "Step",       # Title
    "Thirty", "Five",             # Page number 35
    "Pi",                         # Notation next to page number
    "rules", "and", "lines",      # Collaged words
    "an", "intricate", "spider's", "web", "work"
]

print("--- Syllable Analysis of the Poem ---")
parts, total = count_syllables(poem_words)

print("\n--- Final Calculation ---")
final_equation = " + ".join(parts)
print(f"The total number of syllables is: {final_equation} = {total}")

if total == 17:
    print("\nA poem of 17 syllables is known as an American Sentence.")
else:
    print(f"\nThe poem has {total} syllables, which does not match the American Sentence form.")
