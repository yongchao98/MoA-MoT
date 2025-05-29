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

    def neighbors(word):
        for i in range(len(word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != word[i]:
                    new_word = word[:i] + c + word[i+1:]
                    if is_valid_word(new_word):
                        yield new_word

    queue = deque([(start, [start])])
    visited = set([start])

    while queue:
        current_word, path = queue.popleft()
        if current_word == end:
            return path
        for neighbor in neighbors(current_word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))

    return None

# Find the transformation path
path = transform_word_ladder('REIN', 'CHIS')
if path:
    print(','.join(path))
else:
    print("No valid transformation found.")