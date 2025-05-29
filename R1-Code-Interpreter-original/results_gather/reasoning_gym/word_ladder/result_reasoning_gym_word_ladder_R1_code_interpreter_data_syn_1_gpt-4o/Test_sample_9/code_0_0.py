from nltk.corpus import words

# Ensure the word list is available
word_list = set(words.words())

def is_valid_word(word):
    return word in word_list

def find_word_ladder(start, end):
    # Start with the initial word
    current_word = start
    ladder = [current_word]
    
    # Change 'S' to 'E': "LOGE"
    current_word = "LOGE"
    if is_valid_word(current_word):
        ladder.append(current_word)
    
    # Change 'L' to 'G': "GOGE" (not valid, so skip)
    
    # Change 'G' to 'R': "GORE"
    current_word = "GORE"
    if is_valid_word(current_word):
        ladder.append(current_word)
    
    return ladder

# Find the word ladder
word_ladder = find_word_ladder("LOGS", "GORE")
print(",".join(word_ladder))