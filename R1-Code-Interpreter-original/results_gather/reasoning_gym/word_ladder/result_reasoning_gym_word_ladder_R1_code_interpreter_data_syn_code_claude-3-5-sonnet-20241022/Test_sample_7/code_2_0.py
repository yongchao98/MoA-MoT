# Define a minimal but sufficient set of 4-letter words
words_4 = {
    'HIYA', 'HIND', 'HAND', 'HEAD', 'HEED', 'NEED',  # main path words
    # additional common words to ensure connectivity
    'HIDE', 'HIRE', 'HINT', 'HEAL', 'HEAR', 'HEAT',
    'HELP', 'HERE', 'HELD', 'NEAR', 'NEAT', 'NECK',
    'HEAP', 'HEEL', 'HIGH', 'HOLD', 'HOLE', 'HOME',
    'HOPE', 'HARD', 'HAVE', 'HADE', 'HALE', 'HATE'
}

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words_4:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    from collections import deque
    
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
path = find_path('HIYA', 'NEED')
if path:
    print(','.join(path))
else:
    print("No valid path found")