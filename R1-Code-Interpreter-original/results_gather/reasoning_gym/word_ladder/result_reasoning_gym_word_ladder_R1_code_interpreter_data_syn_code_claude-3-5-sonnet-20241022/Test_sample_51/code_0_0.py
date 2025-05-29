from collections import deque, defaultdict

# First, let's create a set of common 4-letter words
words = set('''
REIN RAIN RAIL RAID PAID PAIN GAIN GRIN CHIN THIN THAN THAT THIS CHIS
MAIN MAIL RAIL FAIL FALL FILL FILM FISH HEAD HEAL HEAR HEIR HAIR PAIR
RAIN REIN REIGN SIGN SINE SITE SITS FITS FINS PINS PINS PIES LIES LIEN
REAL SEAL SEAT BEAT BELT BENT RENT RANT RANK TANK TALK TALL TELL CELL
RITE RATE LATE LAKE LIKE LIME TIME TIME TIRE WIRE WISE RISE RISK RINK
'''.split())

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in words or end not in words:
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

# Find the path
path = find_path('REIN', 'CHIS')
if path:
    print(','.join(path))
else:
    print("No valid path found")