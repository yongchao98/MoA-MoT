from collections import deque, defaultdict

# Define a set of common 4-letter words
WORDS = {
    'PLOY', 'PLAY', 'PLEA', 'PLEA', 'POLE', 'MOLE', 'MOPE',
    'PLOT', 'PLOW', 'PLAN', 'PLEA', 'PLOD', 'PROP', 'PLOP',
    'MOLE', 'POLE', 'POPE', 'PORE', 'POSE', 'POKE', 'PONE',
    'MOPE', 'MOVE', 'MODE', 'MORE', 'MOLE', 'MOTE', 'MOLE'
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
    if start == end:
        return [start]
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                
                if neighbor == end:
                    return new_path
                    
                queue.append(new_path)
    
    return None

# Find the path
path = find_path('PLOY', 'MOPE')
if path:
    print(','.join(path))
else:
    print("No path found")