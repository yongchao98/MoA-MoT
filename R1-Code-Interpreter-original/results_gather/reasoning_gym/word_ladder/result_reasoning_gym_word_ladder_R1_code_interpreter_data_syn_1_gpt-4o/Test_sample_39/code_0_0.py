from nltk.corpus import words

# Ensure the word list is available
try:
    word_list = set(words.words())
except LookupError:
    import nltk
    nltk.download('words')
    word_list = set(words.words())

def is_valid_word(word):
    return word in word_list

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
    
    return []

# Find the transformation path
path = transform_word_ladder('SINS', 'WEST')
print(','.join(path))