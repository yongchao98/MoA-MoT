import itertools

# Load a list of valid English words
with open('/usr/share/dict/words', 'r') as f:
    valid_words = set(word.strip().lower() for word in f)

# Scrambled sentence
scrambled_sentence = "eysa uoY mya ues shit eBook fro nrleay ayn puorpse"

# Function to unscramble a word
def unscramble_word(scrambled_word):
    # Generate all permutations of the scrambled word
    permutations = itertools.permutations(scrambled_word)
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        possible_word = ''.join(perm)
        if possible_word in valid_words:
            return possible_word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word.lower()) for word in scrambled_sentence.split()]

# Capitalize the first word and join the sentence
unscrambled_sentence = ' '.join(unscrambled_words).capitalize()

print(unscrambled_sentence)