from collections import deque

# Expanded set of common 4-letter English words
WORDS = set([
    'DICE', 'DINE', 'FINE', 'FIRE', 'FARE', 'FACE', 'FADE', 'FAME', 'FATE',
    'TEEN', 'TEEM', 'TEAM', 'TEAR', 'TEAL', 'TEAK', 'TEAS', 'TENS', 'TEND',
    'RICE', 'RACE', 'PACE', 'PARE', 'PINE', 'DACE', 'DAME', 'DOME', 'DONE',
    'TONE', 'TUNE', 'TINE', 'TERN', 'TURN', 'TORN', 'TOWN', 'TOON', 'DIES',
    'TIES', 'TEES', 'TIEN', 'TICK', 'TICK', 'TIDE', 'TILE', 'TIME', 'TIRE',
    'VICE', 'VILE', 'VINE', 'WINE', 'WIRE', 'WISE', 'WIFE', 'WIDE', 'TIDE',
    'TILE', 'TINE', 'DIME', 'DIVE', 'DIRE', 'TIRE', 'TIED', 'TIER', 'TEER'
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