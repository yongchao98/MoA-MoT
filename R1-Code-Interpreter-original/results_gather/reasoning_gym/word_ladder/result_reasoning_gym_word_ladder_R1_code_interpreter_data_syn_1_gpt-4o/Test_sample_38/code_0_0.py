from nltk.corpus import words

# Ensure the words corpus is available
import nltk
nltk.download('words')

# Set of valid English words
valid_words = set(words.words())

def is_valid_word(word):
    return word in valid_words

def transform_word_ladder(start, end):
    from collections import deque

    # Queue for BFS
    queue = deque([(start, [start])])
    visited = set([start])

    while queue:
        current_word, path = queue.popleft()

        if current_word == end:
            return path

        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != current_word[i]:
                    next_word = current_word[:i] + c + current_word[i+1:]
                    if next_word not in visited and is_valid_word(next_word):
                        visited.add(next_word)
                        queue.append((next_word, path + [next_word]))

    return None

# Find the transformation path
path = transform_word_ladder('JUGS', 'REFS')
print(','.join(path))