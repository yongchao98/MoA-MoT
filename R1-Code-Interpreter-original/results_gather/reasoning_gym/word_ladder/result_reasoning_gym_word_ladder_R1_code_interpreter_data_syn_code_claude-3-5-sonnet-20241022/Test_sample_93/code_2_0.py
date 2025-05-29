from collections import deque

# Define a carefully curated set of 4-letter English words
# Including more possible intermediate words
WORDS = set([
    'GIVE', 'GOVE', 'GAVE', 'GAPE', 'GATE', 'GORE', 'GONE',
    'LOVE', 'LIVE', 'LIFE', 'LIKE', 'LINE',
    'SAVE', 'SAME', 'SALE', 'SAFE',
    'VASE', 'VALE', 'VOTE', 'VICE',
    'VEEP', 'DEEP', 'DEED', 'VEER',
    'VIVE', 'FIVE', 'DIVE', 'VIDE',
    'VINE', 'VANE', 'VAPE', 'VAMP',
    'CAVE', 'DAVE', 'WAVE', 'RAVE',
    'HAVE', 'HIVE', 'HIDE', 'HIRE',
    'MOVE', 'MORE', 'MODE', 'MADE',
    'PIPE', 'POPE', 'PORE', 'PARE'
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