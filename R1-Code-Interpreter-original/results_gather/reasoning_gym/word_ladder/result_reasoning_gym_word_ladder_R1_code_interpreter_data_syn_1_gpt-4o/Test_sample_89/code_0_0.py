from nltk.corpus import words

# Ensure the word list is available
try:
    word_list = set(words.words())
except LookupError:
    import nltk
    nltk.download('words')
    word_list = set(words.words())

def is_valid_word(word):
    return word.lower() in word_list

def find_word_ladder(start, end):
    from collections import deque

    # Initialize the queue with the start word
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
                    if is_valid_word(next_word) and next_word not in visited:
                        visited.add(next_word)
                        queue.append((next_word, path + [next_word]))

    return None

# Find the word ladder
ladder = find_word_ladder('GHAT', 'TUCK')
print(','.join(ladder) if ladder else "No valid ladder found")