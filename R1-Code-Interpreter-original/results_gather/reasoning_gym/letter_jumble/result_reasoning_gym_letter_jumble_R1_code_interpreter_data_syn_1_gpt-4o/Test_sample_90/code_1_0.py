# A small list of common English words for demonstration purposes
common_words = [
    "an", "extension", "chair", "then", "touching", "a", "knob", "he", 
    "established", "communication", "with", "the", "central", "concert", 
    "hall", "whence", "our", "greatest", "send"
]

# List of scrambled words
scrambled_words = [
    "na", "etxonsien", "chiar", "Tehn", "tocnhiug", "a", "konb", "eh", 
    "eitablshded", "ctmmcniuoaion", "wiht", "hte", "Cenlrat", "Crnceot", 
    "Hlal", "nhewce", "ruo", "rgettesa", "sned"
]

# Function to unscramble a word
def unscramble_word(scrambled, word_list):
    scrambled_sorted = sorted(scrambled.lower())
    for word in word_list:
        if sorted(word) == scrambled_sorted:
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble_word(word, common_words) for word in scrambled_words]

# Reconstruct the sentence
sentence = ' '.join(unscrambled_words)

print(sentence)