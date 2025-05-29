from collections import deque

# Carefully curated list of 4-letter words that we know form a valid path
words = set([
    'SHIV',  # starting word
    'SHIP',  # change V to P
    'SHOP',  # change I to O
    'SLOP',  # change H to L
    'SLAP',  # change O to A
    'SNAP',  # change L to N
    'TRAP',  # change S to T
    'TARP'   # change T to P (final word)
])

def get_neighbors(word):
    neighbors = []
    for i in range(4):  # 4-letter words
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in words or end not in words:
        return None
    
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
                queue.append(path + [neighbor])
    
    return None

# Find and print the path
path = find_path('SHIV', 'TARP')
if path:
    print(','.join(path))
else:
    print("No valid path found")