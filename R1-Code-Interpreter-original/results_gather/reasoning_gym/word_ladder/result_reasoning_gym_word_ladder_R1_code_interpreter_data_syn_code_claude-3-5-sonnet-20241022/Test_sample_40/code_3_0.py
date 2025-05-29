from collections import deque

# Carefully selected word list focusing on common 4-letter words
words = {
    'USED', 'USES', 'USER', 'FUSE', 'FUEL', 'FIRE', 'FILE',
    'URED', 'USE', 'FUME', 'FAME', 'FACE', 'FADE', 'FATE',
    'FUSE', 'FUME', 'FUEL', 'FULL', 'FILL', 'FILE',
    'MUSE', 'MULE', 'MILE', 'MIRE', 'MISE', 'MINE',
    'LOSE', 'LOVE', 'LIVE', 'LIFE', 'LINE', 'LIKE',
    'POSE', 'POLE', 'PILE', 'PINE', 'PIKE', 'PIPE',
    'ROSE', 'ROLE', 'RULE', 'RUSE', 'RILE', 'RIFE'
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