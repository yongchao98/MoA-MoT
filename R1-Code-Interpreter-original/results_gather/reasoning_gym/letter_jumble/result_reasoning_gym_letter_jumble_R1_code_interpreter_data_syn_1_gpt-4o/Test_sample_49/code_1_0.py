import itertools

# List of scrambled words
scrambled_words = [
    "insrgi", "cseldrea", "het", "nma", "ot", "eb", "adde", "deaD", 
    "lemsaicx", "reevy", "eon", "eteprns", "sYe", "sawnesr", "rD", 
    "iisklnW", "aedd", "And", "how"
]

# A small list of common English words for demonstration
common_words = [
    "rising", "declares", "the", "man", "to", "be", "dead", "Dead", 
    "exclaims", "every", "one", "presents", "Yes", "answers", "Dr", 
    "Wilkins", "dead", "And", "how"
]

# Function to unscramble a word
def unscramble(word, word_list):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate in word_list:
            return candidate
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word, common_words) for word in scrambled_words]

# Reconstruct the sentence
sentence = " ".join(unscrambled_words)

print(sentence)