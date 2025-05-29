from collections import deque

# A more comprehensive set of common 4-letter English words
WORDS = {
    'GOLF', 'REFS', 'GOOF', 'GOES', 'GORE', 'GORY', 'GOUT', 'GOAT', 'GOAD',
    'GOLD', 'FOLD', 'FOOD', 'FOOT', 'FORT', 'FORE', 'FIRE', 'FARE', 'FAME',
    'FATE', 'FETE', 'FEDS', 'FEES', 'REFS', 'ROLL', 'ROLE', 'RULE', 'RILE',
    'RIFE', 'RIFT', 'REEF', 'REEK', 'SEEK', 'SEEM', 'STEM', 'STEP', 'STOP',
    'SHOP', 'SHOT', 'SHOE', 'RODE', 'ROME', 'ROPE', 'RIPE', 'RIFE', 'SAFE',
    'SELF', 'SERF', 'SURF', 'TURF', 'ROLF', 'WOLF', 'HALF', 'HOLE', 'HOPE',
    'COPE', 'CORE', 'CURE', 'PURE', 'PORE', 'PORT', 'PART', 'PARK', 'PERK',
    'PEAK', 'PEAL', 'REAL', 'SEAL', 'SEAM', 'TEAM', 'TEAR', 'REAR', 'REEF',
    'GELS', 'GETS', 'LETS', 'LEES', 'FEES', 'FOES', 'TOES', 'TIES', 'LIES',
    'LIED', 'LIED', 'LIDS', 'KIDS', 'KISS', 'MISS', 'MOSS', 'MESS', 'LESS',
    'LENS', 'TENS', 'TENS', 'DENS', 'DENT', 'RENT', 'REST', 'PEST', 'PELT',
    'PELT', 'BELT', 'BOLT', 'BOAT', 'BEAT', 'BEAR', 'REAR', 'REAP', 'LEAP',
    'LEAD', 'LOAD', 'ROAD', 'ROAR', 'REAR', 'FEAR', 'FEAT', 'FEAT', 'FEET',
    'FEED', 'REED', 'REEF', 'BEEF', 'BEER', 'BEAR', 'HEAR', 'HEIR', 'HAIR',
    'FAIR', 'FAIL', 'FALL', 'FILL', 'FILE', 'FIRE', 'FARE', 'FORE', 'FORK',
    'WORK', 'WORD', 'WARD', 'WARN', 'WORN', 'WORM', 'FORM', 'FOAM', 'ROAM',
    'ROAD', 'READ', 'REED', 'FEED', 'FEEL', 'PEEL', 'PEER', 'PEAR', 'FEAR',
    'FEAT', 'FLAT', 'FLAY', 'PLAY', 'PRAY', 'PREY', 'PREP', 'PEEP', 'REEP',
    'RESP', 'REST', 'PEST', 'POST', 'PORT', 'FORT', 'FORD', 'FOLD', 'HOLD',
    'HOLE', 'HOPE', 'ROPE', 'RIPE', 'RIFE', 'LIFE', 'LIFT', 'LEFT', 'LENT',
    'LENS', 'LEES', 'LETS', 'PETS', 'PENS', 'PENT', 'PERT', 'PERL', 'PERK',
    'PEAK', 'PEAR', 'REAR', 'REAP', 'RESP', 'REDS', 'REFS'
}

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

# Find the path
result = find_path('GOLF', 'REFS')
if result:
    print(','.join(result))
else:
    print("No valid path found")