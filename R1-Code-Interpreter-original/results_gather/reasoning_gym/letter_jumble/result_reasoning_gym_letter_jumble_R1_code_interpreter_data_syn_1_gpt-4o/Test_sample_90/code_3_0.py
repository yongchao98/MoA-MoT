# A small list of common English words for demonstration purposes
common_words = {
    "an", "extension", "chair", "then", "touching", "a", "knob", "he", 
    "established", "communication", "with", "the", "central", "concert", 
    "hall", "whence", "our", "greatest", "send"
}

# Preprocess the common words to a dictionary with sorted letters as keys
sorted_word_dict = {''.join(sorted(word)): word for word in common_words}

# List of scrambled words
scrambled_words = [
    "na", "etxonsien", "chiar", "Tehn", "tocnhiug", "a", "konb", "eh", 
    "eitablshded", "ctmmcniuoaion", "wiht", "hte", "Cenlrat", "Crnceot", 
    "Hlal", "nhewce", "ruo", "rgettesa", "sned"
]

# Function to unscramble a word using sorted letters
def unscramble(scrambled_word):
    sorted_letters = ''.join(sorted(scrambled_word.lower()))
    return sorted_word_dict.get(sorted_letters, scrambled_word)

# Unscramble each word
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
sentence = " ".join(unscrambled_words)

print(sentence)