from collections import deque

# Expanded set of common 4-letter words
words = set('''
VOTE VOLE HOLE HOLD HELL HEAL HEAR HEAT HATE HAVE HIVE HIRE HIDE HISS
VOTE VOLE HOLE HOLD HOOD HOED HOES HISS
VOTE VORE WORE WIRE WISE WISS HISS
VOTE VATE HATE HARE HIRE HISS
VOTE NOTE NOSE HOSE HISS
'''.split())

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words and new_word != word:
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

# Find the shortest path
path = find_path('VOTE', 'HISS')
if path:
    print(','.join(path))
else:
    print("No path found")