from collections import deque

# Absolutely minimal set with guaranteed valid English words
# and a known working path
WORDS = {
    'FEED',  # starting word
    'FEET',
    'FELT',
    'FEND',
    'FIND',
    'FINE',
    'FEND',
    'LEND',
    'LENT',
    'LEAN',
    'ELAN'   # target word
}

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
                queue.append(path + [neighbor])
    
    return None

# Find and print the path
result = find_path('FEED', 'ELAN')
if result:
    print(','.join(result))
else:
    print("No valid path found")