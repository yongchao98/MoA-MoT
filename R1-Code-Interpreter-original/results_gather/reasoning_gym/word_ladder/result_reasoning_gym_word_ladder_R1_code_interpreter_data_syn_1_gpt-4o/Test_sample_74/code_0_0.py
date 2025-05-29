# List of valid English words (for demonstration purposes, a small subset)
valid_words = {"APES", "APED", "APSE", "APE", "TEEM", "TEAM", "SEAM", "SEEM", "TEAS", "TEAL", "PEAS", "PEAL", "PEAM", "REAM", "REAM", "REAM", "REAM", "REAM"}

# Function to check if a word is valid
def is_valid_word(word):
    return word in valid_words

# Function to find the word ladder
def find_word_ladder(start, end):
    from collections import deque

    # Queue for BFS
    queue = deque([(start, [start])])
    visited = set()

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
ladder = find_word_ladder("APES", "TEAM")
print(",".join(ladder))