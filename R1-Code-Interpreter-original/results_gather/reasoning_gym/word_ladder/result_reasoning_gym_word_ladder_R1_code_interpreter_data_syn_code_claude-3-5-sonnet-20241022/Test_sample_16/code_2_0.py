from collections import deque

# Further expanded set of 4-letter English words, focusing on possible transition words
WORDS = set([
    'DICE', 'DACE', 'FACE', 'FACT', 'FADE', 'FAME', 'FATE', 'FETE',
    'TEEN', 'SEEN', 'SEED', 'FEED', 'FEET', 'FEEL', 'FELL', 'FILL',
    'RICE', 'RACE', 'PACE', 'PARE', 'PINE', 'DINE', 'FINE', 'FIRE',
    'DIRE', 'TIRE', 'TIER', 'TEER', 'TIEN', 'TEEN',
    'DICK', 'TICK', 'TICK', 'TIDE', 'TILE', 'TIME', 'TINE',
    'DIES', 'TIES', 'TEES', 'TENS', 'DENS', 'DENT', 'TENT',
    'DICE', 'DINE', 'TINE', 'TEEN'
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