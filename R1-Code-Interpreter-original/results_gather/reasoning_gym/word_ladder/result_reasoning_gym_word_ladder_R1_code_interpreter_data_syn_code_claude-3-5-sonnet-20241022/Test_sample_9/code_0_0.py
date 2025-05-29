from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = set([
    'LOGS', 'LOGE', 'LOSE', 'LOVE', 'MOVE', 'MORE', 'GORE',
    'LAGS', 'LEGS', 'LENS', 'LONE', 'LONG', 'LOAD', 'GOAD',
    'GOES', 'GORES', 'LORE', 'WORE', 'WAGE', 'GAGE', 'GAME',
    'COME', 'DOME', 'HOME', 'HOPE', 'HOLE', 'ROLE', 'ROPE',
    'CORE', 'FORE', 'PORE', 'POLE', 'PAGE', 'RAGE', 'RACE'
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

# Find the path and print result
path = find_path('LOGS', 'GORE')
if path:
    print(','.join(path))
else:
    print("No valid path found")