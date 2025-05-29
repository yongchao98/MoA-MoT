from collections import deque

# Create a minimal set of 4-letter words that should work
words = set([
    'JIGS',  # start
    'JAGS',  # change I to A
    'JETS',  # change A to E, G to T
    'JEST'   # change T to S
])

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words:
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

# Find the path
path = find_path('JIGS', 'JEST')
if path:
    print(','.join(path))
else:
    print("No solution found")