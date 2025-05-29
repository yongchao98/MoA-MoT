# List of valid English words for the purpose of this problem
valid_words = set([
    "SCAB", "SLAB", "SLAD", "SLED", "SIDE"
])

# Function to check if two words differ by exactly one letter
def one_letter_diff(word1, word2):
    return sum(a != b for a, b in zip(word1, word2)) == 1

# Starting and ending words
start_word = "SCAB"
end_word = "SIDE"

# Initialize the queue with the start word
queue = [(start_word, [start_word])]

# Perform a breadth-first search
while queue:
    current_word, path = queue.pop(0)
    if current_word == end_word:
        print(",".join(path))
        break
    for word in valid_words:
        if one_letter_diff(current_word, word) and word not in path:
            queue.append((word, path + [word]))