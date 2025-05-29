from itertools import permutations

# List of scrambled words
scrambled_words = ["wlil", "eb", "eknild", "ot", "eht", "rPojtce", "Guerngetb", "Lsceine", "fro", "lal", "osrkw"]

# Manually defined dictionary of possible words
possible_words = [
    "will", "be", "linked", "to", "the", "Project", "Gutenberg", "License", "for", "all", "works"
]

# Function to unscramble a word
def unscramble(scrambled_word, dictionary):
    # Generate all permutations of the scrambled word
    permutations_set = set(permutations(scrambled_word))
    # Check each permutation to see if it's a valid word in the dictionary
    for perm in permutations_set:
        candidate = ''.join(perm)
        if candidate.lower() in map(str.lower, dictionary):
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word, possible_words) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)