from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = set([
    'WAYS', 'WAGS', 'WADS', 'WADS', 'WANS', 'WAPS', 'WARD', 'WARE', 'WARS', 'WARY',
    'WAST', 'WATS', 'WAYS', 'CAYS', 'CODS', 'COGS', 'COKE', 'COLD', 'COLE', 'COLS',
    'COLT', 'COMA', 'COMB', 'COME', 'COMP', 'CONS', 'CONY', 'COPE', 'COPS', 'COPY',
    'CORD', 'CORE', 'CORK', 'CORN', 'CORS', 'CORY', 'COSH', 'COST', 'COSY', 'COTS',
    'WOES', 'GOES', 'DOES', 'COES', 'COLS', 'COYS', 'TOYS', 'TAYS', 'PAYS', 'COTS',
    'WATS', 'CATS', 'COTS', 'COPS', 'COPS', 'CAPS', 'CAYS'
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
path = find_path('WAYS', 'COSY')
if path:
    print(','.join(path))
else:
    print("No path found")