# Load a list of valid English words
with open('/usr/share/dict/words', 'r') as f:
    valid_words = set(word.strip().lower() for word in f)

# Scrambled sentence
scrambled_sentence = "eysa uoY mya ues shit eBook fro nrleay ayn puorpse"

# Function to unscramble a word
def unscramble_word(scrambled_word):
    # Sort the characters of the scrambled word
    sorted_scrambled = ''.join(sorted(scrambled_word))
    # Check each word in the dictionary
    for word in valid_words:
        if ''.join(sorted(word)) == sorted_scrambled:
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word.lower()) for word in scrambled_sentence.split()]

# Capitalize the first word and join the sentence
unscrambled_sentence = ' '.join(unscrambled_words).capitalize()

print(unscrambled_sentence)