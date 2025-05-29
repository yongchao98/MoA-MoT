# Using a very specific set of common 4-letter English words
# These words are all verified to be common English words
WORDS = {
    'PLOY',  # play or scheme (starting word)
    'PLAY',  # to engage in activity for enjoyment
    'PLAN',  # a detailed proposal
    'PLEA',  # an earnest request
    'PLOT',  # a storyline or scheme
    'PLOW',  # to turn up the earth
    'PLOD',  # to move heavily
    'PLUG',  # to fill a hole
    'PLUS',  # positive
    'PREY',  # an animal hunted for food
    'PRAY',  # to communicate with divine power
    'PROP',  # to support
    'PROM',  # a formal dance
    'POLE',  # a long, slender piece of wood
    'POPE',  # head of Catholic Church
    'PORE',  # a minute opening
    'POSE',  # to present oneself
    'POKE',  # to prod
    'MOPE',  # to be dejected (target word)
    'MOLE',  # a small burrowing animal
    'MOVE',  # to change position
    'MODE',  # a way or manner
    'MORE',  # a greater amount
    'HOPE',  # to wish for
    'COPE',  # to deal with
    'ROPE',  # a thick cord
    'DOPE',  # information
}

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS and new_word != word:
                # Check that only one letter is different
                diff = sum(1 for a, b in zip(word, new_word) if a != b)
                if diff == 1:
                    neighbors.append(new_word)
    return neighbors

def find_path(start, end, max_depth=6):
    if start == end:
        return [start]
    
    visited = {start}
    queue = [(start, [start])]
    
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

# Try to find the path
result = find_path('PLOY', 'MOPE')
if result:
    print(','.join(result))
else:
    print("Attempting alternative paths...")
    # Try some known intermediate steps
    for intermediate1 in ['PLAY', 'PLEA', 'PLOT', 'PROP', 'POLE']:
        path1 = find_path('PLOY', intermediate1)
        if path1:
            for intermediate2 in ['POPE', 'HOPE', 'COPE', 'MOLE']:
                path2 = find_path(intermediate1, intermediate2)
                if path2:
                    path3 = find_path(intermediate2, 'MOPE')
                    if path3:
                        # Combine paths, removing duplicates
                        full_path = path1 + path2[1:] + path3[1:]
                        print(','.join(full_path))
                        break