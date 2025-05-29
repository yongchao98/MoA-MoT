from collections import deque

# Expanded set of common 4-letter English words
words = set([
    'FOAL', 'FOOL', 'FOUL', 'FALL', 'FAIL', 'FOLD', 'FOOD', 'FOOT', 'FORD',
    'DROP', 'DRIP', 'PROP', 'CROP', 'DRAP', 'TRAP', 'TRIP', 'TROP', 'FROG',
    'COAL', 'COAT', 'GOAT', 'GOAL', 'ROAD', 'ROAM', 'FOAM', 'FROM', 'FORM',
    'FRAP', 'FLAP', 'FLIP', 'FLOP', 'FLOG', 'FROG', 'PROG', 'PROD', 'PROP',
    'DORP', 'DORP', 'DOLL', 'ROLL', 'ROLE', 'ROPE', 'RIPE', 'FORE', 'FORK',
    'FORT', 'PORT', 'PART', 'DART', 'DIRT', 'DORM', 'DOOM', 'ROOM', 'ROMP',
    'DUMP', 'DAMP', 'DAME', 'DOME', 'DOPE', 'DOER', 'DOES', 'DOSE', 'DOTE',
    'FOAM', 'ROAM', 'ROAD', 'LOAD', 'LOAF', 'LOAN', 'LEAN', 'LEAD', 'FOLD',
    'HOLD', 'HOLE', 'HOPE', 'DOPE', 'ROPE', 'RIPE', 'RAPE', 'RAPS', 'RAMP',
    'ROMP', 'ROME', 'DOME', 'DORM', 'FORM', 'FORT', 'PORT', 'PORE', 'POLE',
    'POLL', 'POOL', 'COOL', 'FOOL', 'FOOD', 'FOLD', 'FOND', 'FONT', 'FORT',
    'FORA', 'FORD', 'FORE', 'FORK', 'FORM', 'FOAM', 'FLAM', 'FLAP', 'FLAT',
    'FIAT', 'FLOAT', 'GLOAT', 'GROAT', 'GROAN', 'GROWN', 'DROWN', 'DRONE',
    'DROVE', 'DRAPE', 'DRAP', 'DROP'
])

def differ_by_one(word1, word2):
    diff = 0
    for c1, c2 in zip(word1, word2):
        if c1 != c2:
            diff += 1
        if diff > 1:
            return False
    return diff == 1

def find_path(start, end):
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        for next_word in words:
            if next_word not in visited and differ_by_one(word, next_word):
                visited.add(next_word)
                queue.append((next_word, path + [next_word]))
    
    return None

# Find the path
path = find_path('FOAL', 'DROP')
if path:
    print(','.join(path))
else:
    print("No path found")