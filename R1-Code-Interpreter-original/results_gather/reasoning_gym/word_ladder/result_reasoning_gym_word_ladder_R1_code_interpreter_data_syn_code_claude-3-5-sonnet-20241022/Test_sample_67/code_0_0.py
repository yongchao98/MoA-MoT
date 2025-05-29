from collections import deque

def get_valid_words():
    # Common 4-letter English words
    words = '''NAGS NOGS FOGS FOPS NAPS NAPE NAPE NAGS NOGS FOGS FOPS
               NABS NOBS FOBS NAPS NAGS NOGS FOGS FOPS
               NAGS NOGS FOGS FOPS'''.split()
    return set(words)

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

def find_path(start, end, valid_words):
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor in valid_words and neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    return None

# Find the path
valid_words = get_valid_words()
path = find_path('NAGS', 'FOPS', valid_words)

if path:
    print(','.join(path))
else:
    print("No valid path found")