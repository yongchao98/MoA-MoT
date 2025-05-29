from collections import deque

# Enhanced 4-letter word list including more common English words
words = set('''
SLOG SLAG SLAT SLOT PLOT PLOD PROD PROM PRIM TRIM TRAM TEAM TEAL TELL
TILL TILL TILT TIDY SLIT SLIP SLOP STOP STEP STEM THEM THEY THEY THEY
TIED TIDY SLAG SLAP SNAP SOAP SOAK SOCK SOCK DOCK DECK DESK DUSK DISK
DINK DINE TINE TIDE TIDY SLIT SHIT SHOT SHOP STOP STEP STEEP STEER STERN
TERN TEEN THEN THIN TIDE TIDY
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

# Find the path from SLOG to TIDY
path = find_path('SLOG', 'TIDY')
if path:
    print(','.join(path))
else:
    print("No valid path found")