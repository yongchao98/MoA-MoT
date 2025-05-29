from collections import deque

# A simple set of valid words for demonstration purposes
valid_words = set([
    "DAZE", "DARE", "DIRE", "DINE", "DINS", "INNS"
])

def is_valid_word(word):
    return word in valid_words

def word_ladder(start, end):
    queue = deque([(start, [start])])
    visited = set([start])
    while queue:
        current_word, path = queue.popleft()
        if current_word == end:
            return path
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                next_word = current_word[:i] + c + current_word[i+1:]
                if next_word != current_word and is_valid_word(next_word) and next_word not in visited:
                    queue.append((next_word, path + [next_word]))
                    visited.add(next_word)  # Avoid revisiting
    return []

# Find the word ladder
ladder = word_ladder("DAZE", "INNS")
print(",".join(ladder))