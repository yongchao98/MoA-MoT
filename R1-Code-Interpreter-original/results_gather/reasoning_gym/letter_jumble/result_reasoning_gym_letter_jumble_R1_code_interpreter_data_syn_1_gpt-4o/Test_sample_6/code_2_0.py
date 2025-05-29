import itertools

# A small sample dictionary of common words for demonstration purposes
common_words = {
    "but", "of", "yesterday", "a", "valuable", "invention", "indeed", "and", "mr", "smith", "this", "morning", "was", "not", "dragging"
}

# Function to unscramble a single word
def unscramble_word(scrambled):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm).lower()
        if candidate in common_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Scrambled sentence
scrambled_sentence = "ubt fo seytyedar A vulalbae enointivn ineedd dna rM miSth shit nnomirg wsa nto nggiard fo"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)