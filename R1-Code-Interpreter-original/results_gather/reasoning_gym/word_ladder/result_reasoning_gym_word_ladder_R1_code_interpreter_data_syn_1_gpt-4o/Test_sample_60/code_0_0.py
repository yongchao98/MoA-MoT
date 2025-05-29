from nltk.corpus import words

# Ensure the word list is available
import nltk
nltk.download('words')

# Set of valid English words
word_set = set(words.words())

def is_valid_word(word):
    return word in word_set

def transform_word_ladder(start, end):
    from collections import deque
    
    # Queue for BFS
    queue = deque([(start, [start])])
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == end:
            return path
        
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != current_word[i]:
                    new_word = current_word[:i] + c + current_word[i+1:]
                    if is_valid_word(new_word) and new_word not in path:
                        queue.append((new_word, path + [new_word]))
    
    return []

# Find the transformation path
path = transform_word_ladder('SLOG', 'TIDY')
print(','.join(path))