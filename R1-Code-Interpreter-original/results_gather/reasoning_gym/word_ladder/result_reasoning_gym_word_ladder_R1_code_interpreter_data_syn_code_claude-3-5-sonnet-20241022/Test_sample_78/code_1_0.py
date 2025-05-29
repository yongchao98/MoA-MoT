from collections import deque

# Define a more comprehensive set of common 4-letter words
WORDS = {
    'VEND', 'VENT', 'VEST', 'VAST', 'HAST', 'HAND', 'HEAD', 'HEAT',
    'WENT', 'WANT', 'WAST', 'BEND', 'BAND', 'SEND', 'LEND', 'LENT',
    'LAST', 'FAST', 'EAST', 'HARD', 'HAVE', 'HATE', 'HALE', 'HALT',
    'VENT', 'VEST', 'BEST', 'BENT', 'BEAT', 'HEAT', 'HEAR', 'HERD',
    'HELD', 'HELP', 'HEAP', 'HEAR', 'HEAD', 'HEAL', 'HAIL', 'HAIR',
    'VELD', 'WELD', 'WEND', 'WIND', 'WING', 'WINE', 'MINE', 'MINT',
    'HINT', 'HUNT', 'HURT', 'HARD', 'CARD', 'CART', 'CAST', 'PAST',
    'POST', 'PORT', 'PART', 'HART', 'HARE', 'HATE', 'HATS', 'HAST'
}

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in WORDS or end not in WORDS:
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
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Find the path
path = find_path('VEND', 'HAST')
if path:
    print(','.join(path))
else:
    print("No path found")