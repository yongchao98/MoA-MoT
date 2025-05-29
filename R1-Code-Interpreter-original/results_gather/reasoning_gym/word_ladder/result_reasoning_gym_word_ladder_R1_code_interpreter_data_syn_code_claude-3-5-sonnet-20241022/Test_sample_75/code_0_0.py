from collections import deque, defaultdict

# First, let's create a set of common 4-letter words
words = set('''
WERE WIRE WIPE WIDE WADE WARE WORE WOVE WAVE WANE WINE WING
TING TRIG TRIP TRAP TRAM TEAM TEAR TEAR TIER TIRE TINE TUNE
TONE TORE WORE WARE WERE WIRE WIPE RIPE RAPE RARE WARE WERE
TRIG TRIP GRIP GRIT GRID GRIM TRIM TRAM TEAM TEAR WEAR WERE
'''.split())

def get_neighbors(word):
    neighbors = []
    for i in range(4):
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

# Find the path from 'WERE' to 'TRIG'
path = find_path('WERE', 'TRIG')
if path:
    print(','.join(path))
else:
    print("No path found")