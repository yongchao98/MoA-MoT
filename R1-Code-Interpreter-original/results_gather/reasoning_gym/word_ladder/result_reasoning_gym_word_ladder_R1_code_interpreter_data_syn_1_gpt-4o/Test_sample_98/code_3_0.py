# Expanded list of valid English words for the purpose of this problem
valid_words = set([
    "HYPE", "TYPE", "TAPE", "TAKE", "TACK", "PACK", "PUCK",
    "HYKE", "HIKE", "LIKE", "LUKE", "LUCK", "PUCK"
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

        # If we reached the end word, return the path
        if current_word == end:
            return path

        # Try changing each letter
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != current_word[i]:
                    next_word = current_word[:i] + c + current_word[i+1:]
                    if is_valid_word(next_word) and next_word not in path:
                        queue.append((next_word, path + [next_word]))

    return []

# Find the word ladder from HYPE to PUCK
ladder = find_word_ladder("HYPE", "PUCK")
print(",".join(ladder))