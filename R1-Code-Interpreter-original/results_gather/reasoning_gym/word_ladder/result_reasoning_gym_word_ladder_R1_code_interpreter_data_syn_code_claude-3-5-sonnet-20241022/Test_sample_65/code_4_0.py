# Using a carefully curated set of common 4-letter English words
WORDS = {
    'PLOY',  # starting word
    'PLAY',  # to engage in activity
    'PLAY',  # variant of play
    'SLAY',  # to kill
    'STAY',  # to remain
    'SLAP',  # to strike
    'SOAP',  # cleaning agent
    'SOUP',  # liquid food
    'SOUR',  # having an acid taste
    'POUR',  # to flow
    'POUT',  # to protrude the lips
    'PORT',  # harbor
    'PART',  # portion
    'PARE',  # to cut off outer covering
    'PORE',  # small opening
    'POLE',  # long rounded piece
    'POPE',  # head of church
    'MOPE',  # target word
    'MOVE',  # to change position
    'MODE',  # manner
    'MOLE',  # small animal
    'MORE',  # greater amount
    'HOPE',  # to wish for
    'COPE',  # to deal with
    'DOPE'   # information
}

def hamming_distance(word1, word2):
    return sum(1 for a, b in zip(word1, word2) if a != b)

def get_neighbors(word):
    return [w for w in WORDS if hamming_distance(word, w) == 1]

def find_shortest_path(start, end, max_depth=10):
    if start == end:
        return [start]
    
    queue = [(start, [start])]
    visited = {start}
    
    while queue:
        current, path = queue.pop(0)
        
        if len(path) > max_depth:
            continue
            
        for next_word in get_neighbors(current):
            if next_word not in visited:
                new_path = path + [next_word]
                
                if next_word == end:
                    return new_path
                    
                visited.add(next_word)
                queue.append((next_word, new_path))
    
    return None

# Find all possible paths
path = find_shortest_path('PLOY', 'MOPE')
if path:
    print(','.join(path))
else:
    # Try specific known good intermediate paths
    for intermediate in ['PLAY', 'PORE', 'POLE', 'POPE']:
        path1 = find_shortest_path('PLOY', intermediate, 3)
        if path1:
            path2 = find_shortest_path(intermediate, 'MOPE', 3)
            if path2:
                full_path = path1 + path2[1:]
                print(','.join(full_path))
                break