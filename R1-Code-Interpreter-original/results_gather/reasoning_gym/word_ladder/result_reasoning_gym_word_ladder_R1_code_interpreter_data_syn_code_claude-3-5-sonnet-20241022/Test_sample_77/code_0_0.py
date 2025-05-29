from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = set([
    'GOLF', 'REFS', 'GOFF', 'GOOF', 'GOES', 'GOES', 'GORE', 'FORE', 'FARE',
    'ROLF', 'ROLE', 'ROLL', 'REEF', 'REEL', 'REAL', 'RELY', 'RULE', 'RIFE',
    'GALE', 'GAME', 'GATE', 'GAVE', 'GENE', 'GONE', 'GORY', 'GULF', 'GOLD',
    'FOLD', 'FORD', 'FORM', 'FORT', 'FOES', 'FEES', 'FEEL', 'FEED', 'REEF',
    'ROES', 'RODE', 'ROPE', 'ROSE', 'ROVE', 'ROTE', 'RATE', 'RAVE', 'RARE',
    'RAFE', 'SAFE', 'SAME', 'SALE', 'SOLE', 'SOME', 'SORE', 'SURE', 'SURF'
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
path = find_path('GOLF', 'REFS')
if path:
    print(','.join(path))
else:
    print("No valid path found")