# Comprehensive list of common 4-letter English words
# Focusing specifically on words that could form a path between PLOY and MOPE
WORDS = {
    'PLOY',  # Starting word
    'PLAY',
    'PLOT',
    'PLOD',
    'PLOP',  # Key intermediate word
    'POOP',
    'POPE',  # Key intermediate word
    'PORE',
    'POLE',
    'MOPE',  # Target word
    'POLY',
    'MOLY',
    'MOLE',
    'MOVE',
    'MODE'
}

def get_neighbors(word):
    neighbors = []
    for i in range(4):  # 4-letter words
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, target):
    if start == target:
        return [start]
    
    # Queue of (word, path)
    queue = [(start, [start])]
    seen = {start}
    
    while queue:
        current, path = queue.pop(0)
        for next_word in get_neighbors(current):
            if next_word not in seen:
                new_path = path + [next_word]
                if next_word == target:
                    return new_path
                seen.add(next_word)
                queue.append((next_word, new_path))
    return None

# Find the path
path = find_path('PLOY', 'MOPE')
if path:
    print(','.join(path))
else:
    print("No valid path found")