from itertools import permutations

def unscramble_word(scrambled, word_list):
    for perm in permutations(scrambled):
        possible_word = ''.join(perm)
        if possible_word.lower() in word_list:
            return possible_word
    return scrambled

# A simple word list for demonstration purposes
word_list = {"your", "discovery", "the"}

scrambled_sentence = "uyor eoicdvrsy Teh"
scrambled_words = scrambled_sentence.split()

unscrambled_words = [unscramble_word(word, word_list) for word in scrambled_words]
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)