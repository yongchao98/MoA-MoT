from collections import deque, defaultdict

# Define a set of common 4-letter words
WORDS = set([
    'GEAR', 'GEAR', 'BEAR', 'BEAT', 'BEAM', 'BEAK', 'BECK', 'DECK', 'DEAR', 'FEAR',
    'HEAR', 'HEAT', 'MEAT', 'MEAL', 'MEAN', 'LEAN', 'LEAP', 'LEAD', 'READ', 'ROAD',
    'ROAM', 'FOAM', 'TEAM', 'TEAR', 'TEAK', 'TEAL', 'TELL', 'TALL', 'TAIL', 'MAIL',
    'MAIN', 'PAIN', 'PAIR', 'PARK', 'PEAK', 'PEAR', 'PEER', 'PEEL', 'POOL', 'COOL',
    'COAL', 'COAT', 'CAME', 'CAMP', 'CAMP', 'CARP', 'CARE', 'CAKE', 'CAVE', 'CAGE',
    'GAME', 'GATE', 'GAZE', 'HAZE', 'MAZE', 'MALE', 'MILE', 'MINE', 'MINT', 'HINT',
    'HEAD', 'HEAL', 'HEAP', 'HELP', 'HELD', 'HOLD', 'HOLE', 'HOME', 'HOPE', 'COPE',
    'CORE', 'CARE', 'CART', 'CAST', 'CASE', 'CASH', 'CALM', 'CAMO'
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

# Find the path from GEAR to CAMO
path = find_path('GEAR', 'CAMO')
if path:
    print(','.join(path))
else:
    print("No valid path found")