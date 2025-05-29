from collections import deque, defaultdict

# Common 4-letter English words
words = '''SHIP SHOP SLIP SLAP SNAP TRAP TARP TRIP TRIM SLIM SLAM SHAM SHAW SHOW SHOT SPOT SPAN SHIN THIN THAN THAT TRAP TARP TEMP TIME TAME TAPE TALE TALK WALK WALL TALL TAIL MAIL MALL MALE MAKE TAKE TAME SAME SAVE HAVE HATE HARE HIRE HIDE SIDE SITE SITS BITS BITE BILE MILE MINE MINT HINT HUNT HURT HURL HULL HALL HALT SALT SALE SAFE SAVE SAKE TAKE TAPE TAPS TOPS TIPS TINS SINS SIPS SITS FITS FILM FILL FALL FAIL TAIL MAIL MALL MALE MAKE TAKE TAME TAPE TARP SHIP SHOP SHOT SHOE SHOD SHED SEED SEND SAND BAND BANK BARK PARK PART TARP TRAP TRIP TRIM GRIM GRID GRIN GAIN PAIN PAIR PARK PERK PEAK PEAR TEAR TSAR TARP TEMP HEMP HELP HELD HOLD HOLE HOPE HOSE ROSE RISE RISK DISK DUSK DESK DECK DOCK ROCK RACK RANK RANT RANT CANT CAST COST COAT CHAT THAT THAN THIN SHIN SHIP SLIP SLAP SNAP SPAN SPAT SPOT SHOT SHOE SHIV'''

# Create a set of valid words
valid_words = set(words.split())

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in valid_words:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in valid_words or end not in valid_words:
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

# Find the path from SHIV to TARP
path = find_path('SHIV', 'TARP')
if path:
    print(','.join(path))
else:
    print("No valid path found")