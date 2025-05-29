from collections import deque

# Expanded word list focusing on possible intermediate words
words = {
    'USED', 'USE', 'FUSE', 'FADE', 'FILE', 'FINE', 'FIRE', 'FIDE',
    'FUEL', 'FEED', 'FLED', 'FOLD', 'FIND', 'FILL', 'FILD',
    'URED', 'USES', 'USER', 'FUED', 'FUND', 'FEND', 'FEED',
    'FLED', 'FREE', 'TREE', 'FLEE', 'FEEL', 'FELL', 'FAIL',
    'FALL', 'FULL', 'FURL', 'FUEL', 'FUME', 'FAME', 'FACE',
    'FADE', 'FATE', 'FETE', 'FILE'
}

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

path = find_path('USED', 'FILE')
if path:
    print(','.join(path))
else:
    print("No solution found")