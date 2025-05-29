# A small list of common English words for demonstration purposes
common_words = {
    "by", "completion", "time", "you", "mean", "yes", "good", "then", 
    "for", "a", "little", "while", "said", "mr", "smith", "as", "he", "deserves"
}

# List of scrambled words
scrambled_words = [
    "yb", "iosnelptCro", "emti", "ouy", "neam", "sYe", "oodG", "yb", "tnhe", 
    "rfo", "a", "tetlli", "hielw", "idas", "rM", "mtSih", "sa", "eh", "ervdees"
]

# Function to unscramble a word
def unscramble_word(scrambled, dictionary):
    from itertools import permutations
    for perm in permutations(scrambled):
        candidate = ''.join(perm).lower()
        if candidate in dictionary:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble_word(word, common_words) for word in scrambled_words]

# Reconstruct the sentence
sentence = ' '.join(unscrambled_words)

print(sentence)