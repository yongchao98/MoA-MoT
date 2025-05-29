from collections import deque

# Load a more comprehensive set of 4-letter words
WORDS = set([
    'GOLF', 'REFS', 'GOOF', 'GOES', 'GORE', 'FORE', 'FARE', 'ROLF', 'ROLE',
    'ROLL', 'REEF', 'REEL', 'REAL', 'RELY', 'RULE', 'RIFE', 'GALE', 'GAME',
    'GATE', 'GAVE', 'GENE', 'GONE', 'GORY', 'GULF', 'GOLD', 'FOLD', 'FORD',
    'FORM', 'FORT', 'FOES', 'FEES', 'FEEL', 'FEED', 'REEF', 'ROES', 'RODE',
    'ROPE', 'ROSE', 'ROVE', 'ROTE', 'RATE', 'RAVE', 'RARE', 'RAFE', 'SAFE',
    'SAME', 'SALE', 'SOLE', 'SOME', 'SORE', 'SURE', 'SURF', 'GOLF', 'GOOF',
    'ROOF', 'ROOK', 'ROCK', 'RACK', 'RACE', 'FACE', 'FAME', 'FATE', 'FETE',
    'FEDS', 'FEES', 'REFS', 'REPS', 'REST', 'RUST', 'RUSE', 'RISE', 'RISK',
    'GELD', 'HELD', 'HELP', 'HEAP', 'HEAD', 'HEAR', 'HEIR', 'HALF', 'HALL',
    'HOLD', 'HOLE', 'HOME', 'HOPE', 'ROPE', 'RIPE', 'RIFE', 'RIFT', 'GIFT'
])

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return None

# Try to find the path
result = find_path('GOLF', 'REFS')
if result:
    print(','.join(result))
else:
    print("No valid path found")