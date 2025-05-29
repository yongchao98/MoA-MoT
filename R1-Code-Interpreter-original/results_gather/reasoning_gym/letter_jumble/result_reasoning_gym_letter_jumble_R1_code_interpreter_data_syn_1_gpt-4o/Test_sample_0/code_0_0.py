from itertools import permutations

# List of known planet names
planet_names = ["Mercury", "Venus", "Mars"]

# Function to find unscrambled word
def unscramble(scrambled_word, word_list):
    scrambled_word_sorted = sorted(scrambled_word.lower())
    for word in word_list:
        if sorted(word.lower()) == scrambled_word_sorted:
            return word
    return None

# Scrambled words
scrambled_words = ["mofr", "uecyrrM", "eunVs", "raMs"]

# Unscramble each word
unscrambled_words = [unscramble(word, planet_names) for word in scrambled_words]

# Construct the sentence
sentence = f"{unscrambled_words[0]} {unscrambled_words[1]} {unscrambled_words[2]} and {unscrambled_words[3]}"
print(sentence)