import itertools

# A small list of common English words for demonstration purposes
word_list = {
    "united", "states", "without", "permission", "and", "paying"
}

# Scrambled sentence
scrambled_sentence = "Utiedn tteaSs wtthoiu erpmisinos dna iwtohut ngiypa"

# Function to unscramble a word
def unscramble_word(scrambled_word):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word.lower()))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate in word_list:
            # Return the word with the original capitalization
            if scrambled_word.istitle():
                return candidate.capitalize()
            else:
                return candidate
    return scrambled_word  # Return the original if no match is found

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)