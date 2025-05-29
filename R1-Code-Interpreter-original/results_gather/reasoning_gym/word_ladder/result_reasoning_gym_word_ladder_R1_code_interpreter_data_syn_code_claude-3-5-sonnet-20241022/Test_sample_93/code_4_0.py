from collections import deque

# Carefully selected word list focusing on very common 4-letter words
WORDS = set([
    'GIVE', 'GAVE', 'GATE', 'GAPE', 'GAME',
    'HAVE', 'HIVE', 'HIDE', 'HIRE',
    'LIVE', 'LOVE', 'LOSE', 'LINE',
    'MOVE', 'MORE', 'MAKE', 'MALE',
    'SAVE', 'SAME', 'SOME', 'SOLE',
    'VALE', 'VOTE', 'VICE', 'VILE',
    'VINE', 'VANE', 'VAPE', 'VEER',
    'VEEP', 'DEEP', 'KEEP', 'PEEP',
    'DIVE', 'FIVE', 'PIPE', 'POPE',
    'RAVE', 'RATE', 'RIPE', 'ROPE'
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

# Find the path and print it
path = find_path('GIVE', 'VEEP')
if path:
    print(','.join(path))
else:
    print("No valid path found")