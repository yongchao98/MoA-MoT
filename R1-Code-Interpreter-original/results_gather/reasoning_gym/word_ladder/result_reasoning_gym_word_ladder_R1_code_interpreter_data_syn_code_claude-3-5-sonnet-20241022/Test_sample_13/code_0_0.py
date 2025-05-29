from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = {
    'PAYS', 'PASS', 'PANS', 'PALS', 'PADS', 'PAYS', 'PAST', 'PATS',
    'PASH', 'PATH', 'SASH', 'SASS', 'SAYS', 'SANS', 'SAPS', 'SAGS',
    'SALE', 'SAME', 'SAVE', 'SAFE', 'SAGE', 'SAKE', 'SALT', 'SAND',
    'PANS', 'PADS', 'PASS', 'PAST', 'PACT', 'PACK', 'PACE', 'PALE',
    'PALM', 'PALL', 'PAIR', 'PAIN', 'PAID', 'PAGE'
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
path = find_path('PAYS', 'SASH')
if path:
    print(','.join(path))
else:
    print("No valid path found")