# Importing necessary libraries
from collections import deque

# A small set of English words for demonstration purposes
# In a real scenario, this should be a comprehensive dictionary of English words
english_words = set([
    "FAVS", "FADS", "SADS", "SUDS", "FUSS", "FUSS", "SUSS", "SUDS"
])

# Function to find the word ladder
def find_word_ladder(start, end, word_set):
    # Queue for BFS
    queue = deque([(start, [start])])
    
    while queue:
        current_word, path = queue.popleft()
        
        # If we reached the end word, return the path
        if current_word == end:
            return path
        
        # Try changing each letter in the current word
        for i in range(len(current_word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if c != current_word[i]:
                    next_word = current_word[:i] + c + current_word[i+1:]
                    if next_word in word_set and next_word not in path:
                        queue.append((next_word, path + [next_word]))
    
    return None

# Find the word ladder
ladder = find_word_ladder("FAVS", "SUDS", english_words)
print(ladder)