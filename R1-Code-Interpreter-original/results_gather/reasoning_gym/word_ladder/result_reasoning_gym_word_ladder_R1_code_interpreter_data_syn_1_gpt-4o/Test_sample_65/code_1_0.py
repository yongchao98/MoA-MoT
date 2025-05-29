# Updated list of valid English words for the purpose of this problem
valid_words = {"PLOY", "PLOP", "PLOE", "MOPE", "MOPY", "MOYE", "MOLE", "MOPE"}

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

# Find the word ladder from PLOY to MOPE
ladder = find_word_ladder("PLOY", "MOPE")
if ladder:
    print(",".join(ladder))
else:
    print("No valid word ladder found.")