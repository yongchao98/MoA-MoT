from itertools import permutations

def unscramble(scrambled_word, possible_words):
    for word in possible_words:
        if sorted(scrambled_word) == sorted(word):
            return word
    return scrambled_word

scrambled_sentence = "bwleo F F Pcotejr"
possible_words = ["below", "elbow", "Project", "F"]

unscrambled_words = []
for word in scrambled_sentence.split():
    unscrambled_word = unscramble(word, possible_words)
    unscrambled_words.append(unscrambled_word)

unscrambled_sentence = ' '.join(unscrambled_words)
print(unscrambled_sentence)