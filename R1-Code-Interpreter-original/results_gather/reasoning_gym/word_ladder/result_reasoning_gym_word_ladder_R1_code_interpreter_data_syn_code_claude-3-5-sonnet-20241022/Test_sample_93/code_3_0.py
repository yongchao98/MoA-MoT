from collections import deque

# Expanded word list with more possible connecting words
WORDS = set([
    'GIVE', 'GAVE', 'GAZE', 'GATE', 'GAPE', 'GAME',
    'HAVE', 'HIVE', 'HIDE', 'HIRE', 'HERE', 'HEAP',
    'LIVE', 'LOVE', 'LOSE', 'LONE', 'LINE', 'LIFE',
    'MOVE', 'MORE', 'MOPE', 'MARE', 'MAKE', 'MALE',
    'SAVE', 'SAME', 'SOME', 'SOLE', 'SALE', 'SAFE',
    'VASE', 'VALE', 'VOTE', 'VOID', 'VICE', 'VILE',
    'VINE', 'VANE', 'VAPE', 'VALE', 'VEER', 'PEER',
    'VEEP', 'DEEP', 'DEED', 'DEER', 'PEEP', 'KEEP',
    'CAVE', 'DAVE', 'DIVE', 'FIVE', 'PIPE', 'POPE',
    'RAVE', 'RATE', 'RIPE', 'ROPE', 'HOPE', 'HEAP',
    'GEVE', 'GENE', 'GONE', 'GORE', 'GERE', 'MERE',
    'VEGE', 'VELE', 'VERE', 'VERY', 'VARY', 'VAMP'
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