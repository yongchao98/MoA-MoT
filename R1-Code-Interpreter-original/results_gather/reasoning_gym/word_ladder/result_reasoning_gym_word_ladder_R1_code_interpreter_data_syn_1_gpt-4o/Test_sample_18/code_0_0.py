from nltk.corpus import words

# Ensure the words corpus is available
import nltk
nltk.download('words')

# Set of valid English words
valid_words = set(words.words())

def is_valid_word(word):
    return word.lower() in valid_words

# Starting and ending words
start_word = "HOPS"
end_word = "MINI"

# Possible sequence of transformations
word_ladder = ["HOPS"]

# Function to find the next word in the ladder
def find_next_word(current_word, target_word):
    for i in range(len(current_word)):
        for char in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            if char != current_word[i]:
                new_word = current_word[:i] + char + current_word[i+1:]
                if is_valid_word(new_word) and new_word not in word_ladder:
                    return new_word
    return None

# Build the word ladder
current_word = start_word
while current_word != end_word:
    next_word = find_next_word(current_word, end_word)
    if next_word:
        word_ladder.append(next_word)
        current_word = next_word
    else:
        break

print(word_ladder)