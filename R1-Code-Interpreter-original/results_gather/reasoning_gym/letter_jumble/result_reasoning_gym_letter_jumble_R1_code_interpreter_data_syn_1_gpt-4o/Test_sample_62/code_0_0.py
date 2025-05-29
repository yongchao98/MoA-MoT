import itertools
from nltk.corpus import words

# Load the list of English words
word_list = set(words.words())

# Scrambled words
scrambled_sentence = "remvoing na eneirt ciyt lal ta oecn fmor eno plaec ot naother Hsi speaicl orpject hda ot od iwth teh"
scrambled_words = scrambled_sentence.split()

# Function to unscramble a word
def unscramble(scrambled_word):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate.lower() in word_list:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)