from collections import deque, defaultdict

# Define a set of common 4-letter words
WORDS = set([
    'HOPE', 'HOPS', 'HIRE', 'HIDE', 'MINE', 'MINI', 'MINT', 'MIND', 'MEND',
    'HIPE', 'HOPS', 'HIPE', 'MIPS', 'MAPS', 'MOPS', 'MINS', 'MINK', 'MINT',
    'MIST', 'MISS', 'HISS', 'HITS', 'HINT', 'MINT', 'MIND', 'MINE', 'MINI',
    'HIND', 'HAND', 'HANS', 'HANS', 'HANS', 'MINS', 'MINS', 'MINE', 'MINI'
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

# Find the path
path = find_path('HOPS', 'MINI')
if path:
    print(','.join(path))
else:
    print("No path found")