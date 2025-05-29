# Minimal but complete word list containing exactly the words we need
words = [
    'CLAN',
    'CRAN',
    'CRAM',
    'DRAM',
    'DRUM',
    'GRUM',
    'GRIM',
    'GRIT',
    'GRIN',
    'RUIN',
    'RUSH'
]

from collections import deque

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        for neighbor in get_neighbors(word):
            if neighbor == end:
                return path + [end]
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    return None

# Convert words list to set for O(1) lookup
word_set = set(words)

# Find the path
path = find_path('CLAN', 'RUSH')

if path:
    print(','.join(path))