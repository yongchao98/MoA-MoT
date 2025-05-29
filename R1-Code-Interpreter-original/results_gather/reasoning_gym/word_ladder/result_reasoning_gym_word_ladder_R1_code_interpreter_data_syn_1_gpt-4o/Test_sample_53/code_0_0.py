# List of valid English words for the purpose of this problem
valid_words = set([
    "WAYS", "WARS", "WARY", "CARS", "CARY", "COSY", "CORY", "CAYS", "CAYS"
])

# Function to check if a word is valid
def is_valid_word(word):
    return word in valid_words

# Function to find the word ladder
def find_word_ladder(start, end):
    from collections import deque

    # Queue for BFS
    queue = deque([(start, [start])])

    while queue:
        current_word, path = queue.popleft()

        if current_word == end:
            return path

        # Try changing each letter
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != current_word[i]:
                    new_word = current_word[:i] + c + current_word[i+1:]
                    if is_valid_word(new_word) and new_word not in path:
                        queue.append((new_word, path + [new_word]))

    return None

# Find the word ladder from WAYS to COSY
ladder = find_word_ladder("WAYS", "COSY")
print(",".join(ladder))