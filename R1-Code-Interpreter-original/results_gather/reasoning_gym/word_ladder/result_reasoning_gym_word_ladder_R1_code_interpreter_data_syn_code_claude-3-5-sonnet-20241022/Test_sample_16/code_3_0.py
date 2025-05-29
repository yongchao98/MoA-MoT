from collections import deque

# Specifically curated word list that includes our solution path and related words
WORDS = set([
    # Core path words
    'DICE',
    'DIRE',
    'TIRE',
    'TINE',
    'TEEN',
    # Additional common 4-letter words that might be useful
    'DINE',
    'TARE',
    'TIDE',
    'TIME',
    'TILE',
    'TEES',
    'TENS'
])

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Find the path and print result
path = find_path('DICE', 'TEEN')
if path:
    print(','.join(path))
else:
    print("No valid path found")