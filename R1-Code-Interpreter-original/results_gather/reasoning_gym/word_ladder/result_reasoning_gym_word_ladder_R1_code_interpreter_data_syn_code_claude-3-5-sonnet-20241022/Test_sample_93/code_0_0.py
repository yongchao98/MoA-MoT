from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = set([
    'GIVE', 'GOVE', 'GAVE', 'GATE', 'GAPE', 'GALE', 'GAME', 'GORE', 'GONE',
    'LOVE', 'LIVE', 'LIFE', 'LIKE', 'LINE', 'LIME', 'LANE', 'LATE', 'LAKE',
    'MAKE', 'MALE', 'MILE', 'MINE', 'MOVE', 'MODE', 'MORE', 'MOPE', 'MADE',
    'SAVE', 'SAME', 'SALE', 'SAFE', 'SAGE', 'SANE', 'SAKE', 'SIDE', 'SITE',
    'VASE', 'VALE', 'VALE', 'VOTE', 'VOLE', 'VICE', 'VILE', 'VINE', 'VIBE',
    'VEEP', 'DEEP', 'DEED', 'DEER', 'VEER', 'PEER', 'PEEL', 'PEEK', 'PEEP'
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

# Find the path and print it
path = find_path('GIVE', 'VEEP')
if path:
    print(','.join(path))
else:
    print("No valid path found")