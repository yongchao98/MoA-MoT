# A small sample dictionary of common English words
word_list = {
    "removing", "an", "entire", "city", "all", "at", "once", "from", "one", "place", "to", "another", 
    "his", "special", "project", "had", "to", "do", "with", "the"
}

# Scrambled words
scrambled_sentence = "remvoing na eneirt ciyt lal ta oecn fmor eno plaec ot naother Hsi speaicl orpject hda ot od iwth teh"
scrambled_words = scrambled_sentence.split()

# Function to unscramble a word
def unscramble(scrambled_word):
    for word in word_list:
        if sorted(scrambled_word.lower()) == sorted(word):
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)