# Define a set of common 4-letter words
word_set = {
    'VANS', 'VANE', 'VINE', 'VINS', 'VIES', 'VICE', 'VILE', 'VALE', 'VARS',
    'KANS', 'KINS', 'KIDS', 'KITS', 'KIRS', 'KISS', 'KIPS', 'KITE', 'KINE',
    'PANS', 'PINS', 'PIRS', 'PIES', 'PIER', 'PITS', 'PIKE', 'PILE', 'PINE',
    'MANS', 'MINS', 'MIRS', 'MITE', 'MILE', 'MINE', 'MIKE', 'MIRE', 'MANE',
    'RANS', 'RINS', 'RIPS', 'RIPE', 'RITE', 'RISE', 'RISK', 'RIFT', 'RAVE',
    'SANS', 'SINS', 'SIRS', 'SITE', 'SIRE', 'SIDE', 'SINE', 'SAVE', 'SALE'
}

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    from collections import deque
    
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
path = find_path('VANS', 'KIRS')
if path:
    print(','.join(path))
else:
    print("No solution found")